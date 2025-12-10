import requests
from typing import Dict, List, Any, Optional


class EnrichmentResult:
    """
    Wrapper for raw enrichment analysis output.

    This class stores parsed metadata and the list of raw result entries.
    """

    def __init__(self, enrichmentData: Dict[str, Any]) -> None:
        """
        Parse enrichment output dictionary.

        Parameters
        ----------
        enrichmentData : dict
            Raw enrichment JSON/dictionary, typically as returned by an
            external enrichment tool.
        """
        self.data: Dict[str, Any] = enrichmentData

        self.input: Dict[str, Any] = self.data["results"]["input_list"]
        self.organism: str = self.input["organism"]
        self.mapped_ids: str = self.input["mapped_ids"]
        self.mapped_count: int = self.input["mapped_count"]
        self.unmapped_ids: str = self.input["unmapped_ids"]
        self.unmapped_count: int = self.input["unmapped_count"]

        # Sorted list of enrichment results (typically by FDR)
        self.result: List[Dict[str, Any]] = self.data["results"]["result"]


def prepare_list_for_enrichment(nodes: List[str]) -> str:
    """
    Prepare a list of node identifiers for enrichment analysis by converting
    it into a comma-separated string without quotes or brackets.

    This function takes a list like:
        ["A", "B", "C"]
    and converts it into:
        A, B, C

    Parameters
    ----------
    nodes : list[str]
        List of node identifiers (e.g., gene names).

    Returns
    -------
    str
        A formatted string suitable for use in an enrichment API request.
    """
    as_string: str = str(nodes)[1:-1]
    as_string = as_string.replace("'", "")
    return as_string


def prepare_enrichment_result(enrichment: Dict[str, Any]) -> Optional[EnrichmentResult]:
    """
    Prepare and validate an enrichment result before wrapping it into an
    EnrichmentResult object.

    If the enrichment response contains an error in the 'search' field,
    the function returns None.

    Parameters
    ----------
    enrichment : dict
        Raw enrichment result dictionary returned by the API.

    Returns
    -------
    EnrichmentResult or None
        Wrapped enrichment result if valid, otherwise None.
    """
    if (
        isinstance(enrichment, dict)
        and "search" in enrichment
        and isinstance(enrichment["search"], dict)
        and "error" in enrichment["search"]
    ):
        return None

    enrichment_result: EnrichmentResult = EnrichmentResult(enrichment)
    return enrichment_result


def get_enrichment(
    input_genes_string: str,
    organism_id: str,
    goterm_type: str,
    test_type: str = "FISHER",
    correction: str = "FDR",
) -> Optional[Dict[str, Any]]:
    """
    Request functional enrichment data from the PANTHER database API.

    This function performs a GET request to the PANTHER enrichment endpoint
    using the given genes, organism, and GO term category.

    Parameters
    ----------
    input_genes_string : str
        Comma-separated list of genes.
    organism_id : str
        Organism identifier used by the PANTHER database.
    goterm_type : str
        GO term category: 
        - "MF" for Molecular Function
        - "BP" for Biological Process
        - "CC" for Cellular Component
    test_type : str, optional
        Statistical test type (default is "FISHER").
    correction : str, optional
        Multiple testing correction method (default is "FDR").

    Returns
    -------
    dict or None
        Parsed JSON response from the API if successful, otherwise None.
    """
    match goterm_type:
        case "MF":
            data_type: str = "GO:0003674"
        case "BP":
            data_type: str = "GO:0008150"
        case "CC":
            data_type: str = "GO:0005575"
        case _:
            print('Wrong goterm_type. Use "MF","BP","CC" instead. Ending get_enrichment function.')
            return None

    req_link: str = (
        f"https://pantherdb.org/services/oai/pantherdb/enrich/overrep"
        f"?geneInputList={input_genes_string}"
        f"&organism={organism_id}"
        f"&annotDataSet={data_type}"
        f"&enrichmentTestType={test_type}"
        f"&correction={correction}"
    )

    headers: Dict[str, str] = {"Content-Type": "application/json"}
    response = requests.get(req_link, headers=headers)

    if response.status_code == 200:
        data: Dict[str, Any] = response.json()
        return data

    print("Failed to get data. Ending get_enrichment function.")
    return None


def get_evaluated_nodes(phenotype: List[str], evaluation: bool = True) -> List[str]:
    """
    Filter phenotype nodes based on their evaluation sign.

    The function assumes that each node is prefixed with:
        "+" for positive regulation
        "-" for negative regulation

    If evaluation is True, only positive nodes are returned.
    If evaluation is False, only negative nodes are returned.

    Parameters
    ----------
    phenotype : list[str]
        List of phenotype nodes with "+" or "-" prefixes.
    evaluation : bool, optional
        If True, return positively evaluated nodes.
        If False, return negatively evaluated nodes.
        Default is True.

    Returns
    -------
    list[str]
        List of filtered node names without their prefix.
    """
    result_list: List[str] = []

    for node in phenotype:
        if node[0] == "+" and evaluation:
            result_list.append(node[1:])
        elif node[0] == "-" and not evaluation:
            result_list.append(node[1:])

    return result_list
