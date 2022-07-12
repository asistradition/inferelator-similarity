import pandas as pd
from typing import List

def build_joint_sets(
    bd: pd.DataFrame,
    include_species: bool = True,
    add_to_list: List[set] = None
) -> List[set]:
    """
    Build a list of sets, where each set is highly similar genes,
    from a blast dataframe

    :param bd: Preprocessed blast data
    :type bd: pd.DataFrame
    :param add_to_list: A list of sets to add to, defaults to None
    :type add_to_list: list(set()), optional
    :return: A list of sets, where sets are highly similar genes
    :rtype: list(set())
    """

    valid_genes = {}
    result_dict = {}

    if add_to_list is not None:
        for s in add_to_list:
            for g in s:
                result_dict[g] = s
                valid_genes[g] = True

    # Make dict of sets where each set is group of homologs
    # and the elements of the set all point to the same set
    # Stash the sets in a list for easy access later
    for qgene, sgene, qspecies, sspecies in zip(
        bd['qseqid'],
        bd['sseqid'],
        bd['qspecies'],
        bd['sspecies']
    ):

        i, j = (qspecies, qgene), (sspecies, sgene)

        joint = set([i, j])

        try:
            result_dict[i].update(joint)
        except KeyError:
            result_dict[i] = joint

        try:
            result_dict[j].update(joint)
        except KeyError:
            result_dict[j] = joint

        result_dict[i].update(result_dict[j])
        result_dict[j] = result_dict[i]

        # Make sure result_dict for all the elements in this set are fused
        for k in result_dict[i].copy():
            try:
                if id(result_dict[i]) != id(result_dict[k]):
                    result_dict[i].update(result_dict[k])
            except KeyError:
                pass

            result_dict[k] = result_dict[i]

        valid_genes[i] = True
        valid_genes[j] = True

    # One last unification pass
    for i in list(result_dict.keys()):
        for k in result_dict[i].copy():

            try:
                if id(result_dict[i]) != id(result_dict[k]):
                    result_dict[i].update(result_dict[k])
            except KeyError:
                pass

            result_dict[k] = result_dict[i]

    # Make a list of unique sets from the duplicate references
    # Use a dict with ids to keep track of what we've seen before
    result_list = []
    set_list = {}
    for _, x in result_dict.items():
        try:
            if set_list[id(x)]:
                continue
        except KeyError:
            result_list.append(x)
            set_list[id(x)] = True

    if include_species:
        return result_list

    else:
        return [set(x[1] for x in s) for s in result_list]
