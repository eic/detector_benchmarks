import pandas as pd
import argparse
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser(prog='Plot output branch sizes', description='Plot output branch sizes')

parser.add_argument("-c", dest="current_campaign_file", action="store", required=True, help="Enter the current campaign file")
parser.add_argument("-d", dest="default_file", action="store", required=True, help="Enter the default file")
                    
args=parser.parse_args()


campaign1=args.current_campaign_file
campaign2=args.default_file


# Load the data from the CSV file
df1 = pd.read_csv(campaign1+'.txt', header=None)
df2 = pd.read_csv(campaign2+'.txt', header=None)

# Plot the third column ('Value') against the first column ('Object')
plt.figure(figsize=(10,6))
plt.scatter(df1.iloc[:,0], df1.iloc[:,2])
plt.scatter(df2.iloc[:,0], df2.iloc[:,2])

plt.title("Branch Sizes (Bytes) vs Branch Names")




# Show the figure
plt.tight_layout()
plt.yscale('log')
plt.savefig(campaign1+'_vs_'+campaign2+'.png')

print(df1)
print(df2)

# Assuming both dataframes have the same structure and the first column is branch name
# Merge the two dataframes on the branch name (first column)
merged_df = pd.merge(df1.iloc[:, [0, 2]], df2.iloc[:, [0, 2]], on=df1.columns[0], suffixes=('_' + campaign1, '_' + campaign2))

# Create a new column that calculates the difference between the third columns of the two DataFrames
merged_df['Difference'] = merged_df.iloc[:, 1] - merged_df.iloc[:, 2]

# Create a new DataFrame with the branch names and the difference
result_df = merged_df[[df1.columns[0], 'Difference']]

# Display the resulting DataFrame
print(result_df)

# Sort the DataFrame by the absolute value of the difference in descending order
sorted_df = result_df.reindex(result_df['Difference'].abs().sort_values(ascending=False).index)

# Pick the top 10 branches with the largest differences
top_20_branches = sorted_df.head(20)

# Display the top 10 branches
print(top_20_branches)


# Optionally, save it to a new CSV file
sorted_df.to_csv(f"{campaign1}_vs_{campaign2}_difference.csv", index=False)
