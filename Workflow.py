

# from traceback import print_tb
# from typing import Tuple
# from isdaapi import Protein, Structure, Mutation
# import pandas as pd

# with open("idmapping_2025_12_06.list") as L:
#     with open("Protein_with_Stru_details.tsv", "w") as W:
#         count = 0
#         for line in L:
#             if count >= 1: continue
#             print(count)
#             protein_detail = Protein.get_protein_details(line.strip())
#             if protein_detail.get("StructuralRecords") != 0:
#                 print(f"UniProtID\t{protein_detail.get('uniprot_id')}\tStructural Details\t{protein_detail.get('StructuralRecords')}", file = W)
#             count = count + 1 


# ###P07550


# ##################################################
# """
#  Get the Protein Details
# """
# print(Protein.get_protein_details("P07550"))   ### Add the Annotation Score in the API

# """
#  Get the UniProt Protein sequence 
# """
# Seq = Protein.Sequence("P07550")
# print(Seq)
# """
# #  Fetch the Structural Details
# # """
# print(Structure.StructureDetails("P07550"))
# #print(Structure.StructureDetails("P07550", PDBDetails=True))

# """
# Fetch the Natural Mutation information
# """

# mut_df = Mutation.MutationTable("P07550", source="clinvar", record_type="full")
# mut_df.to_csv(f"MutationRecords_P07550.csv", index=False)


# ## Visulize the Mutation Matrix
# #Mutation.plot_mutation_matrix(mut_df)

# ## Get the Basic Mutation Details
# Mut_Summary = Mutation.GetMutationSummary(mut_df, top_n=20)
# print(pd.DataFrame(Mut_Summary['Mutation_type']))
# print(pd.DataFrame(Mut_Summary['significance_status']).T)

# ## Get the Signifiance level 
# print(Mutation.MatchSignificance(mut_df, pattern="Pathogenic", case=True))

# ## Analyze change in Physiochemcial Mutation Properties

# print(Mutation.AnalyzeMutationProperties(mut_df).head(20))

# """
# Analysis the Disease status for this Protein ID
# """
# print(mut_df['Phenotype/Disease'].value_counts(sort=True))


# ## Generate the Mutated Protein Sequence

# Mutation_list = mut_df[["Alt_AA","Position"]]
# Protein.MutateSequence(Seq, Mutation_list)



# #############################################
# """
# SASA Analysis
# """

# # from isdaapi import SASA

# # sasa_result = SASA.calculate_sasa_for_chain("8uo1", "R")
# # if sasa_result:
# #     ids, names, sasa_vals = sasa_result
# #     for J in range(len(ids)):
# #         print(ids[J], names[J], sasa_vals[J])
    



################################################
""""
Interaction Analysis
"""
from isdaapi import InteractionMining

df_interactions = InteractionMining.analyze_residue_interactions(
    structure_path="8uo1.cif",              ### Provide the CIF file
    selection_list=['/A/20/']
)
print(df_interactions.head())


# ##################################################
# """
# Computational Mutations
# """


# df_ComMutation = pd.DataFrame(Mutation.AFMutationTable(uniprot_id="P07550", ranges=None, significance_type=None))

# """
# pathogenic, ambiguous, benign
# """

# df_ComMutation_withRanges = pd.DataFrame(Mutation.AFMutationTable(uniprot_id="P07550", ranges='45:50', significance_type="pathogenic"))

# print(df_ComMutation_withRanges)

# """
# Map the Computational Mutation Score with the Experimentally known position
# """

# # df_ComMutation_singe = pd.DataFrame(Afmissense.AFMutationTable(uniprot_id="P07550", ranges='292', significance_type=None))



