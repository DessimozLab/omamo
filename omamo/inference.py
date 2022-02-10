import collections
import itertools
import logging
import statistics as stats
import numpy
import pandas
from concurrent.futures import ProcessPoolExecutor, as_completed
from pyoma.browser.db import Database, GeneNameOrSymbolIdMapper, BestIdPerEntryOnlyMixin


logger = logging.getLogger(__name__)
OrthologsOverlapBioProcess = collections.namedtuple("OrthologsOverlapBioProcess",
                                                    ["ortholog", "processes", "similarity"])

SUMMARY_DTYPE = [("GOnr", "i4"), ("Species", "S5"), ("NrOrthologs", "i4"),
                 ("FuncSim_Mean", "f8"), ("FuncSim_Std", "f8"), ("Score", "f8")]
DETAIL_DTYPE = [("GOnr", "i4"), ("Species", "S5"), ("Ref", "i1"), ("EntryNr", "i4"),
                ("Label", "S50")]
IC_DTYPE = [("GOnr", "i4"), ("ic", "f8")]


def find_orthologs(db: Database, query_species: str, model_species: str):
    orthologs = db.get_vpairs_between_species_pair(query_species, model_species)
    return orthologs


def go_overlap(db, orthologs, inf_content, inf_threshold=5):
    """This function determines common GO terms.
     Input: orthologs, a numpy array of pairwise orthologs
     Output is a list of lists:
     [[pair, [bio process GO overlap], perc_similarity],...]"""

    go = db.gene_ontology
    BP = 1  # aspect biological_process

    def get_go_anno_with_parents(entry_nr):
        annos = set(a['TermNr'] for a in db.get_gene_ontology_annotations(entry_nr))
        if len(annos) <= 0:
            return annos
        all_annos = set.union(*(go.get_superterms_incl_queryterm(term) for term in annos))
        return all_annos

    output = []
    for pair in orthologs:
        annos1, annos2 = (get_go_anno_with_parents(pair[z]) for z in ('EntryNr1', 'EntryNr2'))
        # overlap and union contain term that have common/union annotations (any ontology) and a
        # information content >= inf_threshold
        overlap = frozenset(t for t in annos1.intersection(annos2)
                            if inf_content.get(t.id, 0) >= inf_threshold)
        union = frozenset(t for t in annos1.union(annos2)
                          if inf_content.get(t.id, 0) >= inf_threshold)

        bp_terms = frozenset(t for t in overlap if t.aspect == BP)
        if len(bp_terms) > 0:
            intrsctn_inf_cont = sum((inf_content.get(t.id, 0) for t in overlap))
            union_inf_cont = sum((inf_content.get(t.id, 0) for t in union))
            try:
                similarity = intrsctn_inf_cont / union_inf_cont
                p = OrthologsOverlapBioProcess(ortholog=(int(pair['EntryNr1']), int(pair['EntryNr2'])),
                                               processes=[t.id for t in bp_terms],
                                               similarity=similarity)
                output.append(p)
            except ZeroDivisionError:
                pass
    return output


def filter_high_inf_content_orthologs(it, threshold):
    """yields those orthologs that have a GO similarity score >= threshold"""
    yield from (ortholog for ortholog in it if ortholog.similarity >= threshold)


def extract_processes_with_occurency_range(similar_orthologs, min_cnt=0, max_cnt=-1):
    """returns all the processes that occur min_cnt <= x <= max_cnt in the orthologs"""
    counts = collections.Counter(itertools.chain.from_iterable(
        (o.processes for o in similar_orthologs)))
    processes = frozenset((proc for proc, cnt in counts.items()
                           if cnt >= min_cnt and max_cnt < 0 or cnt <= max_cnt))
    return processes


class BestGeneId(GeneNameOrSymbolIdMapper, BestIdPerEntryOnlyMixin):
    def map_entry_nr_range(self, start, stop):
        arr = super().map_entry_nr_range(start, stop)
        return self.filter_best_id(arr)


def get_gene_names(db, similar_orthologs):
    def enr_rng(k):
        s = set(z.ortholog[k] for z in similar_orthologs)
        return min(s), max(s) + 1

    query_range, target_range = (enr_rng(k) for k in (0, 1))
    mapper = BestGeneId(db)
    query_xref = mapper.map_entry_nr_range(query_range[0], query_range[1])
    target_xref = mapper.map_entry_nr_range(target_range[0], target_range[1])
    return {int(e['EntryNr']): e['XRefId'] for e in itertools.chain(query_xref, target_xref)}


def pivot_go_process(db, similar_orthologs, model_species):
    summary = []
    detail = []

    go_processes = extract_processes_with_occurency_range(similar_orthologs, max_cnt=5000)
    gene_names = get_gene_names(db, similar_orthologs)
    sp = model_species.encode('utf-8')
    for go in go_processes:
        similarity, seen_genes = [], set([])
        for sim_orth in similar_orthologs:
            if go in sim_orth.processes:
                for ref, gene in enumerate(sim_orth.ortholog):
                    if gene not in seen_genes:
                        detail.append((go, sp, ref, gene, gene_names[gene],))
                        seen_genes.add(gene)
                similarity.append(sim_orth.similarity)
        avg = stats.mean(similarity)
        std = stats.stdev(similarity) if len(similarity) > 1 else 0
        sim = sum(similarity)
        summary.append((go, sp, len(similarity), avg, std, sim))
    summary = numpy.array(summary, dtype=SUMMARY_DTYPE)
    details = numpy.array(detail, dtype=DETAIL_DTYPE)
    return summary, details


def compute_omamo_for_species(db_path, ic, query_species, model_species):
    with Database(db_path) as db:
        orthologs = find_orthologs(db, query_species, model_species)
        similar_orthologs = list(
            filter_high_inf_content_orthologs(go_overlap(db,
                                                         orthologs=orthologs,
                                                         inf_content=ic),
                                              threshold=0.05))
        return pivot_go_process(db, similar_orthologs, model_species)


def combine_datasets(summaries, details):
    summary = pandas.concat([pandas.DataFrame(z) for z in summaries], ignore_index=True)
    detail = pandas.concat([pandas.DataFrame(z) for z in details], ignore_index=True)
    summary = summary.sort_values(by=['GOnr', 'Score'], ignore_index=True, ascending=[True, False])
    detail = detail.sort_values(by=["GOnr", "Species", "Ref"], ignore_index=True, ascending=True)
    return summary, detail


def write_hdf5(fpath, summary: pandas.DataFrame, detail: pandas.DataFrame, ic: dict):
    import tables
    filters = tables.Filters(complevel=6, complib="zlib")
    with tables.open_file(fpath, 'w', filters=filters) as h5:
        sum_data = summary.to_records(index=False, column_dtypes=dict(SUMMARY_DTYPE))
        h5.create_table('/omamo', 'Summary',
                        title="Summary information of best modal organism per GO process",
                        obj=sum_data,
                        expectedrows=len(sum_data),
                        createparents=True)

        detail_data = detail.to_records(index=False, column_dtypes=dict(DETAIL_DTYPE))
        h5.create_table('/omamo', 'detail',
                        title="Ortholog information for GO-term in species",
                        obj=detail_data,
                        expectedrows=len(detail_data))
        ic_data = numpy.fromiter(ic.items(), dtype=IC_DTYPE)
        h5.create_table('/', 'ic', title="Information Content (IC) of GO terms",
                        obj=ic_data, expectedrows=len(ic_data))
    with tables.open_file(fpath, 'a') as h5:
        sum_tab = h5.get_node('/omamo/Summary')
        detail_tab = h5.get_node('/omamo/detail')
        for col in ("GOnr", "Species", "Score"):
            sum_tab.colinstances[col].create_csindex()
        for col in ("GOnr", "Species", "EntryNr", "Label"):
            detail_tab.colinstances[col].create_csindex()


def write_csv(fpath, summary: pandas.DataFrame, detail: pandas.DataFrame):
    summary['Species'] = summary['Species'].apply(bytes.decode)
    detail['Species'] = detail['Species'].apply(bytes.decode)
    detail['Label'] = detail['Label'].apply(bytes.decode)
    genes = detail.groupby(by=["GOnr", "Species", "Ref"])[['Label']].agg(";".join)
    genes = genes.reset_index().pivot(index=["GOnr", "Species"], columns="Ref", values="Label").reset_index()\
        .rename(columns={0: "QuerySpeciesGenes", 1: "ModelSpeciesGenes"})
    result = pandas.merge(summary, genes, on=["GOnr", "Species"])
    result.sort_values(by=["GOnr", "Score"], ignore_index=True, ascending=[True, False])
    result.to_csv(fpath,
                  columns=["GOnr", "Species", "QuerySpeciesGenes", "ModelSpeciesGenes",
                           "NrOrthologs", "FuncSim_Mean", "FuncSim_Std", "Score"],
                  float_format="%.4f",
                  index=False,
                  sep="\t")


def build_omamo(db_path, ic, query, models, h5_out=None, tsv_out=None):
    with ProcessPoolExecutor() as executor:
        futures = [executor.submit(compute_omamo_for_species, db_path, ic, query, model) for model in models]
        summaries = []
        details = []
        for f in as_completed(futures):
            s, d = f.result()
            summaries.append(s)
            details.append(d)
    summary, detail = combine_datasets(summaries, details)
    if h5_out is not None:
        write_hdf5(h5_out, summary, detail, ic)
    if tsv_out is not None:
        write_csv(tsv_out, summary, detail)
