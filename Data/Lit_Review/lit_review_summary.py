import pandas as pd
import matplotlib.pyplot as plt


#======= Tcell =======
# Load the T-cell dataset
tcell_df = pd.read_csv("bcb330_lit_review_tcell_verified.csv")

# Normalize public access info
tcell_df["public"] = tcell_df["public"].str.lower().replace("application", "no")

# Drop duplicates based on Author, Year, and Sequence Type
tcell_unique = tcell_df.drop_duplicates(subset=["first_author", "Year", "data_type"])

# Group by data type and accessibility
tcell_counts = tcell_unique.groupby(["data_type", "public"]).size().unstack(fill_value=0)

# Plot
tcell_counts.plot(kind="bar", stacked=True, figsize=(10, 6))
plt.title("T-cell Data Types by Public Availability")
plt.ylabel("Number of Studies")
plt.xlabel("T-cell Data Type")
plt.xticks(rotation=45)
plt.tight_layout()
plt.savefig("tcell_data_availability.png")
plt.close()


# ======= Tumor =======
# Load the Tumor dataset
tumour_df = pd.read_csv("bcb330_lit_review_tumour_verified.csv")

# Normalize public access info
tumour_df["tumour_public"] = tumour_df["tumour_public"].str.lower().replace("application", "no")

# Drop duplicates based on Author, Year, and Sequence Type
tumour_unique = tumour_df.drop_duplicates(subset=["first_author", "Year", "data"])

# Group by data type and accessibility
tumour_counts = tumour_unique.groupby(["data", "tumour_public"]).size().unstack(fill_value=0)

# Plot
tumour_counts.plot(kind="bar", stacked=True, figsize=(10, 6))
plt.title("Tumor Data Types by Public Availability")
plt.ylabel("Number of Studies")
plt.xlabel("Tumor Data Type")
plt.xticks(rotation=45)
plt.tight_layout()
plt.savefig("tumour_data_availability.png")
plt.close()


#======= Cancer Type =======
# Count by cancer type
tumour_unique = tumour_df.drop_duplicates(subset=["first_author", "Year", "type"])
tumour_cancer_counts = tumour_unique["type"].value_counts()

# Plot
tumour_cancer_counts.plot(kind="bar", figsize=(10, 6))
plt.title("Cancer Types in Tumor Studies")
plt.ylabel("Number of Studies")
plt.xlabel("Cancer Type")
plt.xticks(rotation=45)
plt.tight_layout()
plt.savefig("tumour_cancer_types.png")
plt.close()
