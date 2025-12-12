
from typing import List, Any

from biodivine_aeon import *

from EnrichmentClasses import (
    EnrichmentPSBN,
    EnrichmentPSBNInstance,
    EnrichmentAttractor
)

from Enrichment import (
    get_evaluated_nodes,
    prepare_list_for_enrichment,
    get_enrichment,
    prepare_enrichment_result,
)

from Visualization import (
    append_column_to_xlsx
)


def pipeline(
    psbn: EnrichmentPSBN,
    network: Any,
    reference_genome_id: str,
    net_name: str,
) -> None:
    """
    Run the enrichment analysis pipeline on a logical network.

    For each color configuration of the network, this function:
      - builds the corresponding fully specified network,
      - computes stable phenotypes and attractors,
      - runs GO Biological Process enrichment for each attractor,
      - writes enrichment results to Excel files (per attractor, per instance,
        and across all instances),
      - and stores the attractors in the provided PSBN enrichment object.

    Parameters
    ----------
    psbn : EnrichmentPSBN
        PSBN enrichment structure that collects all instances created
        in this pipeline.
    network : Any
        Asynchronous graph of the network from provided file.
    reference_genome_id : str
        Identifier of the reference genome used for enrichment (e.g. a
        PANTHER / organism ID).
    net_name : str
        Base name used as a prefix for all generated Excel files.

    Returns
    -------
    None
    """
    stg_general = AsynchronousGraph(network)
    all_colors = stg_general.mk_unit_colors()

    for color_index, color in enumerate(all_colors):
        fully_specified_network = color.instantiate(stg_general.reconstruct_network())
        ctx = SymbolicSpaceContext(fully_specified_network)
        stg = AsynchronousGraph(fully_specified_network, ctx)

        classification = Classification.classify_stable_phenotypes(ctx, stg)
        attractors = Attractors.attractors(stg)
        attractorClassifs = Classification.classify_attractor_bifurcation(stg, attractors)
        attractors_types = list(attractorClassifs)[0].feature_list()

        instance_results: List[List[str]] = []
        for res in range(len(classification)):
            a = list(classification)[res]
            instance_results.append(a.feature_list())

        new_psbn_instance = EnrichmentPSBNInstance()
        i: int = 0

        for nodes, attractor_type in zip(instance_results, attractors_types):
            to_enrich = get_evaluated_nodes(nodes)
            to_enrich = prepare_list_for_enrichment(to_enrich)

            enrichment = get_enrichment(to_enrich, reference_genome_id, "BP")
            enrichment_result = prepare_enrichment_result(enrichment)

            calculated_attractor = EnrichmentAttractor(to_enrich, attractor_type, enrichment_result, 0.05)
            new_psbn_instance.add_attractor(calculated_attractor)

            append_column_to_xlsx(f"{net_name}_OnAttractors.xlsx", calculated_attractor.goterms, column_name=(f"[clr:{color_index}][att:{i}]"))
            i += 1

        new_psbn_instance.set_color(color)

        append_column_to_xlsx(f"{net_name}_OnInstance.xlsx", new_psbn_instance.goterm_id_intersection(), column_name=(f"[{color_index}]"))
        psbn.add_instance(new_psbn_instance)

    append_column_to_xlsx(f"{net_name}_OnAllInstances.xlsx", psbn.goterms_id_intersection_on_all_instances(), column_name="[whole]")
