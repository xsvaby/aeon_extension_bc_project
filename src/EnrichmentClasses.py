from __future__ import annotations
from typing import Dict, List, Set, Any, Optional

from Enrichment import (
    EnrichmentResult,
    get_evaluated_nodes,
    prepare_list_for_enrichment,
    get_enrichment,
    prepare_enrichment_result,
)


class EnrichmentPSBN:
    """
    Represents a collection of PSBN enrichment instances.

    This class allows comparison of GO terms (by ID or name) across multiple
    `EnrichmentPSBNInstance` objects and can compute their intersections.
    """

    def __init__(self, instances: List["EnrichmentPSBNInstance"]) -> None:
        """
        Initialize the PSBN with a list of PSBN instances.

        Parameters
        ----------
        instances : list[EnrichmentPSBNInstance]
            List of enrichment PSBN instances to manage.
        """
        self.instances: List[EnrichmentPSBNInstance] = [instance for instance in instances]
        # Dictionary of all GO terms across all instances: GO ID -> EnrichmentGOterm
        self.all_goterms: Dict[str, EnrichmentGOterm] = {}

    def add_instance(self, instance: "EnrichmentPSBNInstance") -> None:
        """
        Add a new PSBN instance and update the global GO term dictionary.

        Parameters
        ----------
        instance : EnrichmentPSBNInstance
            Instance to be added.
        """
        self.instances.append(instance)
        self.all_goterms.update(instance.get_all_goterms())

    def get_all_goterms(self) -> Dict[str, "EnrichmentGOterm"]:
        """
        Get all GO terms collected from all instances.

        Returns
        -------
        dict[str, EnrichmentGOterm]
            Mapping of GO IDs to GO term objects.
        """
        return self.all_goterms

    def goterms_id_intersection_on_all_instances(self) -> Set[str]:
        """
        Compute the intersection of GO term IDs across all PSBN instances.

        Returns
        -------
        set[str]
            Set of GO term IDs that are common to all instances.
        """
        intersected: Set[str] = self.instances[0].goterm_id_intersection()

        for i, instance in enumerate(self.instances):
            if i == 0:
                continue
            intersected = intersected.intersection(instance.goterm_id_intersection())

        return intersected

    def goterms_name_intersection_on_all_instances(self) -> Set[str]:
        """
        Compute the intersection of GO term names across all PSBN instances.

        Returns
        -------
        set[str]
            Set of GO term names that are common to all instances.
        """
        intersected: Set[str] = set(self.instances[0].goterm_name_intersection())

        for i, instance in enumerate(self.instances):
            if i == 0:
                continue
            intersected = intersected.intersection(instance.goterm_name_intersection())

        return intersected

    def goterms_intersection_on_all_instances(self) -> Dict[str, "EnrichmentGOterm"]:
        """
        Get the GO term objects that are common across all PSBN instances.

        Returns
        -------
        dict[str, EnrichmentGOterm]
            Mapping of GO IDs to GO term objects that appear in every instance.
        """
        intersected_ids: Set[str] = self.goterms_id_intersection_on_all_instances()
        intersected_goterms: Dict[str, EnrichmentGOterm] = {}

        for go_id in intersected_ids:
            intersected_goterms[go_id] = self.get_all_goterms()[go_id]

        return intersected_goterms

    def count_go_ids_frequencies_in_all_instances(self) -> Dict[str, int]:
        """
        Count how often each GO term appears across all PSBN instances.

        This aggregates GO-term frequencies by:
        1. Calling `count_go_ids_frequencies()` on each `EnrichmentPSBNInstance`
        2. Summing the frequencies of identical GO terms across all instances

        Returns
        -------
        dict[str, int]
            Mapping of GO ids to their total frequency across all
            instances and their attractors.
        """
        frequencies: Dict[str, int] = {}

        for instance in self.instances:
            instance_frequencies = instance.count_go_ids_frequencies()

            for go_id, count in instance_frequencies.items():
                frequencies[go_id] = frequencies.get(go_id, 0) + count

        return frequencies
    
    def count_goterms_frequencies_in_all_instances(self) -> Dict[EnrichmentGOterm, int]:
        """
        Convert aggregated GO-term ID frequencies across all instances into
        actual GO-term object frequencies.

        This method:
        1. Retrieves GO-term ID frequencies using `count_go_ids_frequencies_in_all_instances()`
        2. Converts each GO ID into its corresponding `EnrichmentGOterm` object
        using `get_all_goterms()`
        3. Returns a dictionary keyed by GO-term objects instead of IDs

        Returns
        -------
        dict[EnrichmentGOterm, int]
            Mapping from GO-term objects to their total frequency across all instances.
        """
        frequencies_go_ids: Dict[str, int] = self.count_go_ids_frequencies_in_all_instances()
        frequencies_goterms: Dict[EnrichmentGOterm, int] = {}

        for go_id, frequency in frequencies_go_ids.items():
            goterm: EnrichmentGOterm = self.get_all_goterms()[go_id]
            frequencies_goterms[goterm] = frequency

        return frequencies_goterms

    def print_unmapped_ids_per_instance_per_attractor(self) -> None:
        """
        Print unmapped IDs for each attractor within each PSBN instance.

        The output is grouped first by instance (with its color), and then
        by attractor index within that instance.
        """
        for i, instance in enumerate(self.instances):
            print(f"{i}: {instance.color}")

            for j, attractor in enumerate(instance.attractors):
                print(f"{attractor.unmapped_ids} [{j}]")
            
            print("---------")
    
    def print_only_output_unmapped_ids_per_instance_per_attractor(self, output_nodes: Set[str]) -> None:
        """
        Print output unmapped IDs for each attractor within each PSBN instance.

        The output is grouped first by instance (with its color), and then
        by attractor index within that instance.
        """
        for i, instance in enumerate(self.instances):
            print(f"{i}: {instance.color}")

            for j, attractor in enumerate(instance.attractors):
                filtered_to_only_outputs = [node for node in attractor.unmapped_ids_set if node in output_nodes]
                print(f"{filtered_to_only_outputs} [{j}]")
            
            print("---------")


    def unmapped_ids_intersection_on_all_instances(self) -> Set[str]:
        """
        Compute the intersection of unmapped IDs across all PSBN instances.

        The intersection is computed by repeatedly intersecting the unmapped
        ID sets from each instance.

        Returns
        -------
        set[str]
            Set of unmapped IDs that appear in all instances.
        """
        intersected: Set[str] = {}

        for i, instance in enumerate(self.instances):
            if i == 0:
                intersected = instance.unmapped_ids_intersection()
                continue

            intersected = intersected.intersection(
                instance.unmapped_ids_intersection()
            )

        return intersected

    def _count_ids_frequencies_in_all_instances(self, method_name: str) -> Dict[str, int]:
        """
        Internal helper to aggregate ID frequencies across all instances.

        Parameters
        ----------
        method_name : str
            Name of the method on each instance that returns a frequency dict.

        Returns
        -------
        dict[str, int]
            Aggregated mapping of IDs to their total frequency across all instances.
        """
        frequencies: Dict[str, int] = {}

        for instance in self.instances:
            instance_frequencies = getattr(instance, method_name)()

            for _id, freq in instance_frequencies.items():
                frequencies[_id] = frequencies.get(_id, 0) + freq

        return frequencies

    def count_unmapped_ids_frequencies_in_all_instances(self) -> Dict[str, int]:
        """
        Count how often each unmapped ID appears across all PSBN instances.

        The frequencies from each instance are aggregated into a single
        global frequency dictionary.

        Returns
        -------
        dict[str, int]
            Mapping of unmapped IDs to their total frequency across all instances.
        """
        return self._count_ids_frequencies_in_all_instances("count_unmapped_ids_frequencies")

    def count_mapped_ids_frequencies_in_all_instances(self) -> Dict[str, int]:
        """
        Count how often each mapped ID appears across all PSBN instances.

        The frequencies from each instance are aggregated into a single
        global frequency dictionary.

        Returns
        -------
        dict[str, int]
            Mapping of mapped IDs to their total frequency across all instances.
        """
        return self._count_ids_frequencies_in_all_instances("count_mapped_ids_frequencies")

    def count_attractors(self) -> int:
        """
        Count the total number of attractors across all PSBN instances.

        Returns
        -------
        int
            Total number of attractors contained in all instances.
        """
        count: int = 0

        for instance in self.instances:
            count += len(instance.attractors)

        return count


class EnrichmentPSBNInstance:
    """
    Represents a single PSBN instance.

    A PSBN instance contains multiple attractors (`EnrichmentAttractor`)
    and their associated GO terms. It can compute intersections and
    uniqueness of GO terms across those attractors.
    """

    def __init__(self) -> None:
        """
        Initialize an empty PSBN instance.
        """
        self.attractors: List[EnrichmentAttractor] = []
        self.attractor_types: List[str] = []
        self.color = None
        # All GO terms across the attractors in this instance: GO ID -> EnrichmentGOterm
        self.all_goterms: Dict[str, EnrichmentGOterm] = {}

    def add_attractor(self, attractor: "EnrichmentAttractor") -> None:
        """
        Add an attractor to the PSBN instance and update GO terms.

        Parameters
        ----------
        attractor : EnrichmentAttractor
            Attractor to be added.
        """
        self.attractors.append(attractor)
        self.attractor_types.append(attractor.attractor_type)
        self.all_goterms.update(attractor.get_all_goterms())

    def set_color(self, color: str) -> None:
        """
        Assign a display color to this instance, usually for plotting.

        Parameters
        ----------
        color : str
            Color value (e.g., hex code or a named color).
        """
        self.color = color

    def get_all_goterms(self) -> Dict[str, "EnrichmentGOterm"]:
        """
        Get all GO terms associated with this instance.

        Returns
        -------
        dict[str, EnrichmentGOterm]
            Mapping of GO IDs to GO term objects.
        """
        return self.all_goterms

    def goterm_id_intersection(self) -> Set[str]:
        """
        Compute the intersection of GO term IDs shared by all attractors.

        Returns
        -------
        set[str]
            Set of GO term IDs common to all attractors in this instance.
        """
        intersect: Set[str] = self.attractors[0].go_terms_set

        for i, attractor in enumerate(self.attractors):
            if i == 0:
                continue
            intersect = intersect.intersection(attractor.go_terms_set)

        return intersect

    def goterm_name_intersection(self) -> Set[str]:
        """
        Compute the intersection of GO term names shared by all attractors.

        Returns
        -------
        set[str]
            Set of GO term names common to all attractors in this instance.
        """
        intersect: Set[str] = self.attractors[0].get_goterm_labels()

        for i, attractor in enumerate(self.attractors):
            if i == 0:
                continue
            intersect = intersect.intersection(attractor.get_goterm_labels())

        return intersect

    def goterm_intersection(self) -> Dict[str, "EnrichmentGOterm"]:
        """
        Get GO term objects that appear in all attractors.

        Returns
        -------
        dict[str, EnrichmentGOterm]
            Mapping of GO IDs to GO term objects that appear in every attractor.
        """
        intersected_ids: Set[str] = self.goterm_id_intersection()
        intersected_goterms: Dict[str, EnrichmentGOterm] = {}

        for go_id in intersected_ids:
            intersected_goterms[go_id] = self.get_all_goterms()[go_id]

        return intersected_goterms

    def count_go_ids_frequencies(self) -> Dict[str, int]:
        """
        Count how often each GO term appears across all attractors.

        Returns
        -------
        dict[str, int]
            Mapping of GO id to the number of attractors that contain it.
        """
        frequencies: Dict[str, int] = {}

        for go_id, _ in self.all_goterms.items():
            for attractor in self.attractors:
                if go_id in attractor.go_terms_set:
                    frequencies[go_id] = frequencies.get(go_id, 0) + 1

        return frequencies
    
    def count_goterms_frequencies(self) -> Dict[EnrichmentGOterm, int]:
        """
        Convert GO-term ID frequencies within a single instance into
        GO-term object frequencies.

        This method:
        1. Retrieves GO-term ID frequencies using `count_go_ids_frequencies()`
        2. Looks up each GO ID in the instance's GO-term dictionary
        (`get_all_goterms()`)
        3. Returns a mapping keyed by GO-term objects rather than their IDs

        Returns
        -------
        dict[EnrichmentGOterm, int]
            Mapping from GO-term objects to their frequency across attractors
            in this instance.
        """
        frequencies_go_ids: Dict[str, int] = self.count_go_ids_frequencies()
        frequencies_goterms: Dict[EnrichmentGOterm, int] = {}

        for go_id, frequency in frequencies_go_ids.items():
            goterm: EnrichmentGOterm = self.get_all_goterms()[go_id]
            frequencies_goterms[goterm] = frequency

        return frequencies_goterms

    def unmapped_ids_intersection(self) -> Set[str]:
        """
        Compute the intersection of unmapped IDs across all attractors
        within a single PSBN instance.

        Returns
        -------
        set[str]
            Set of unmapped IDs that appear in every attractor.
        """
        intersect: Set[str] = {}

        for i, attractor in enumerate(self.attractors):
            if i == 0:
                intersect = attractor.unmapped_ids_set
                continue

            intersect = intersect.intersection(attractor.unmapped_ids_set)

        return intersect

    def _count_id_frequencies(self, attr_name: str) -> Dict[str, int]:
        """
        Internal helper to count frequencies of IDs across attractors.

        Parameters
        ----------
        attr_name : str
            Name of the attribute on each attractor that holds a set of IDs.

        Returns
        -------
        dict[str, int]
            Mapping of IDs to the number of attractors in which they appear.
        """
        frequencies: Dict[str, int] = {}

        for attractor in self.attractors:
            ids_set = getattr(attractor, attr_name)
            for id in ids_set:
                frequencies[id] = frequencies.get(id, 0) + 1

        return frequencies

    def count_unmapped_ids_frequencies(self) -> Dict[str, int]:
        """
        Count how often each unmapped ID appears across attractors
        in a single PSBN instance.

        Each attractor contributes at most one count per unmapped ID.

        Returns
        -------
        dict[str, int]
            Mapping of unmapped IDs to the number of attractors in which they appear.
        """
        return self._count_id_frequencies("unmapped_ids_set")

    def count_mapped_ids_frequencies(self) -> Dict[str, int]:
        """
        Count how often each mapped ID appears across attractors
        in a single PSBN instance.

        Each attractor contributes at most one count per mapped ID.

        Returns
        -------
        dict[str, int]
            Mapping of mapped IDs to the number of attractors in which they appear.
        """
        return self._count_id_frequencies("mapped_ids_set")

    def __str__(self) -> str:
        result: str = "PSBN instance \n"
        for attractor in self.attractors:
            result += " |-- " + str(attractor) + "\n"
        return result

    def __repr__(self) -> str:
        return self.__str__()


class EnrichmentAttractor:
    """
    Represents a single biological attractor and its enrichment results.

    Contains:
    - A set of enriched nodes (genes, etc.)
    - GO terms satisfying an FDR threshold
    - Basic mapping information from the enrichment tool
    """

    def __init__(self, enriched_nodes: str, attractor_type: str,
                 enrichment_result: Optional["EnrichmentResult"],
                 fdr: float,) -> None:
        """
        Initialize an enrichment attractor.

        Parameters
        ----------
        enriched_nodes : str
            Comma-separated string of node names (e.g., gene symbols).
        attractor_type : str
            Label describing the attractor type (e.g. "upregulated", "cluster 1").
        enrichment_result : EnrichmentResult or None
            Result of an enrichment analysis. If None, no GO terms are added.
        fdr : float
            False Discovery Rate threshold for filtering GO terms.
        """
        self.fdr: float = fdr
        self.attractor_type: str = attractor_type

        # GO ID -> EnrichmentGOterm
        self.goterms: Dict[str, EnrichmentGOterm] = {}
        # Set of GO IDs
        self.go_terms_set: Set[str] = set()
        # Cleaned set of enriched node names
        self.enriched_nodes: Set[str] = {name.strip() for name in enriched_nodes.split(",") if name.strip()}

        self.mapped_ids: str = ""
        self.unmapped_ids: str = ""
        self.unmapped_ids_set: Set[str] = set()
        self.mapped_ids_set: Set[str] = set()

        if enrichment_result is None:
            return

        self.mapped_ids = enrichment_result.mapped_ids
        self.unmapped_ids = enrichment_result.unmapped_ids
        self.unmapped_ids_set = {unmapped for unmapped in self.unmapped_ids.split(",")}
        self.mapped_ids_set = {mapped for mapped in self.mapped_ids.split(",")}


        # Populate GO terms from enrichment result, with FDR filter
        for process in enrichment_result.result:
            go_term = EnrichmentGOterm(process)

            # Skip terms above FDR threshold or "invalid" ones (starting with "-")
            if go_term.fdr > self.fdr or go_term.process_name.startswith("-"):
                continue

            self.goterms[go_term.go_id] = go_term
            self.go_terms_set.add(go_term.go_id)

    def get_goterm_labels(self) -> Set[str]:
        """
        Get all GO term names associated with this attractor.

        Returns
        -------
        set[str]
            Set of GO term names.
        """
        return {goterm.process_name for goterm in self.goterms.values()}

    def get_goterms_by_set(self, wanted: Set[str]) -> List["EnrichmentGOterm"]:
        """
        Return GO term objects for a given set of GO IDs.

        Parameters
        ----------
        wanted : set[str]
            Set of GO IDs of interest.

        Returns
        -------
        list[EnrichmentGOterm]
            GO term objects whose IDs are present in `wanted`.
        """
        return [self.goterms[go_id] for go_id in wanted if go_id in self.goterms]

    def get_all_goterms(self) -> Dict[str, "EnrichmentGOterm"]:
        """
        Return all GO terms for this attractor.

        Returns
        -------
        dict[str, EnrichmentGOterm]
            Mapping of GO IDs to GO term objects.
        """
        return self.goterms

    def __str__(self) -> str:
        return f"{self.attractor_type}"

    def __repr__(self) -> str:
        return f"{self.attractor_type}"


class EnrichmentGOterm:
    """
    Represents a single Gene Ontology (GO) term from enrichment analysis.
    """

    def __init__(self, process: Dict[str, Any]) -> None:
        """
        Initialize GO term from a single enrichment result entry.

        Parameters
        ----------
        process : dict
            Dictionary representing one GO term from the enrichment result.
            Expected keys include: "term", "fold_enrichment", "fdr",
            "expected", "number_in_reference", "pValue", "plus_minus".
        """
        self.go_id: str = process.get("term", {}).get("id", "")
        self.process_name: str = process["term"]["label"]
        self.fold_enrichment: float = process["fold_enrichment"]
        self.fdr: float = process["fdr"]
        self.expected: float = process["expected"]
        self.number_in_reference: int = process["number_in_reference"]
        self.p_value: float = process["pValue"]
        self.plus_minus: str = process["plus_minus"]

        # Relationships in the GO graph:
        # children/parents: EnrichmentGOterm -> relationship string
        self.children: Dict["EnrichmentGOterm", str] = {}
        self.parents: Dict["EnrichmentGOterm", str] = {}

    def add_child(self, child: "EnrichmentGOterm", relation: str) -> None:
        """
        Add a child GO term with a relationship.

        Parameters
        ----------
        child : EnrichmentGOterm
            Child GO term.
        relation : str
            Relationship label (e.g. "part_of", "regulates").
        """
        self.children[child] = relation

    def add_parent(self, parent: "EnrichmentGOterm") -> None:
        """
        Add a parent GO term.

        Parameters
        ----------
        parent : EnrichmentGOterm
            Parent GO term.
        """
        self.parents[parent] = ""

    def __repr__(self) -> str:
        return f"{self.plus_minus}{self.process_name}"

    def __str__(self) -> str:
        return f"{self.plus_minus}{self.process_name}"
