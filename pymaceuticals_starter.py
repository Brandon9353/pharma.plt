#!/usr/bin/env python
# coding: utf-8

# # Pymaceuticals Inc.
# ---
# 
# ### Analysis
# 
# - Add your analysis here.
#  

# In[1]:


# Dependencies and Setup
import matplotlib.pyplot as plt
import pandas as pd
import scipy.stats as st

# Study data files
mouse_metadata_path = "data/Mouse_metadata.csv"
study_results_path = "data/Study_results.csv"

# Read the mouse data and the study results
mouse_metadata = pd.read_csv(mouse_metadata_path)
study_results = pd.read_csv(study_results_path)


# In[2]:


mouse_metadata


# In[3]:


study_results


# In[4]:


# Combine the data into a single DataFrame
mouse_study = pd.merge(mouse_metadata, study_results, how="left", on=["Mouse ID", "Mouse ID"])

# Display the data table for preview
mouse_study.head()


# In[5]:


# Checking the number of mice.
len(mouse_study["Mouse ID"].unique())


# In[6]:


mouse_study


# In[7]:


# Our data should be uniquely identified by Mouse ID and Timepoint
# Get the duplicate mice by ID number that shows up for Mouse ID and Timepoint.
duplicates = mouse_study[mouse_study.duplicated(subset=["Mouse ID", "Timepoint"])]["Mouse ID"]

duplicates


# In[8]:


# Optional: Get all the data for the duplicate mouse ID.
duplicates_II = mouse_study[mouse_study["Mouse ID"] == "g989"]
duplicates_II


# In[9]:


# Create a clean DataFrame by dropping the duplicate mouse by its ID.
clean_mouse_study = mouse_study[mouse_study["Mouse ID"].isin(duplicates) == False]
clean_mouse_study


# In[10]:


# Checking the number of mice in the clean DataFrame.
len(clean_mouse_study["Mouse ID"].unique())


# ## Summary Statistics

# In[11]:


summary_stats = clean_mouse_study.groupby("Drug Regimen")["Tumor Volume (mm3)"].agg(['mean', 'median', 'var', 'std', 'sem'])
summary_stats


# In[12]:


# Generate a summary statistics table of mean, median, variance, standard deviation, and SEM of the tumor volume for each regimen

# Use groupby and summary statistical methods to calculate the following properties of each drug regimen:
# mean, median, variance, standard deviation, and SEM of the tumor volume.
summary_stats = clean_mouse_study.groupby("Drug Regimen")["Tumor Volume (mm3)"].agg(["mean", "median", "var", "std", "sem"])


summary_df = pd.DataFrame({
    "Mean": summary_stats["mean"],
    "Median": summary_stats["median"],
    "Variance": summary_stats["var"],
    "Std Deviation": summary_stats["std"],
    "SEM": summary_stats["sem"]
})
summary_df


# In[13]:


# A more advanced method to generate a summary statistics table of mean, median, variance, standard deviation,
# and SEM of the tumor volume for each regimen (only one method is required in the solution)

# Using the aggregation method, produce the same summary statistics in a single line
mean = clean_mouse_study.groupby("Drug Regimen")["Tumor Volume (mm3)"].mean()
median = clean_mouse_study.groupby("Drug Regimen")["Tumor Volume (mm3)"].median()
variance = clean_mouse_study.groupby("Drug Regimen")["Tumor Volume (mm3)"].var()
std = clean_mouse_study.groupby("Drug Regimen")["Tumor Volume (mm3)"].std()
sem = clean_mouse_study.groupby("Drug Regimen")["Tumor Volume (mm3)"].sem()
pd.DataFrame(
    {
      "Mean Tumor Volume" : mean,
      "Median Tumor Volume" : median,
      "Tumor Volume Variance" : variance,
      "Tumor Volume Std. Dev" : std,
      "Tumor Volume Std. Err" : sem  
    }
)


# In[14]:


# Generate a bar plot showing the total number of rows (Mouse ID/Timepoints) for each drug regimen using Pandas.
counts = clean_mouse_study["Drug Regimen"].value_counts()
counts.plot(kind="bar")
plt.xlabel("Drug Regimen")
plt.xticks(rotation = 50)
plt.ylabel("number of mice tested")
plt.show


# In[15]:


# Generate a bar plot showing the total number of rows (Mouse ID/Timepoints) for each drug regimen using pyplot.
counts = clean_mouse_study["Drug Regimen"].value_counts()
plt.bar(counts.index.values,counts.values)
plt.xlabel("Drug Regimen")
plt.xticks(rotation = 50)
plt.ylabel("number of mice tested")
plt.show


# In[16]:


# Generate a pie chart, using Pandas, showing the distribution of unique female versus male mice used in the study

# Get the unique mice with their gender


# Make the pie chart
counts = clean_mouse_study.Sex.value_counts()
counts.plot(kind = "pie", autopct = "%1.1f%%")


# In[17]:


# Generate a pie chart, using pyplot, showing the distribution of unique female versus male mice used in the study

# Get the unique mice with their gender


# Make the pie chart
counts = clean_mouse_study.Sex.value_counts()
plt.pie(counts.values, labels = counts.index.values, autopct = "%1.1f%%")
plt.ylabel("Sex")
plt.show


# ## Quartiles, Outliers and Boxplots

# In[18]:


# Calculate the final tumor volume of each mouse across four of the treatment regimens:
# Capomulin, Ramicane, Infubinol, and Ceftamin

# Start by getting the last (greatest) timepoint for each mouse
max_tumor = clean_mouse_study.groupby(["Mouse ID"])["Timepoint"].max()
max_tumor = max_tumor.reset_index()
# Merge this group df with the original DataFrame to get the tumor volume at the last timepoint
mouse_tumor = max_tumor.merge(clean_mouse_study, on= ["Mouse ID", "Timepoint"], how="left")
mouse_tumor


# In[19]:


# Put treatments into a list for for loop (and later for plot labels)
treatment_list = ["Capomulin", "Ramicane", "Infubinol", "Ceftamin"]

# Create an empty list to hold tumor volume data for all treatments
tumor_vol_list = []

# Iterate through each drug regimen
for drug in treatment_list:
    # Subset the data to get tumor volume for the current drug regimen
    final_tumor_vol = mouse_tumor.loc[mouse_tumor["Drug Regimen"] == drug, "Tumor Volume (mm3)"]
    
    # Calculate quartiles and IQR
    quartiles = final_tumor_vol.quantile([0.25, 0.5, 0.75])
    lowerq = quartiles[0.25]
    upperq = quartiles[0.75]
    iqr = upperq - lowerq
    lower_bound = lowerq - (1.5 * iqr)
    upper_bound = upperq + (1.5 * iqr)
    
    # Filter outliers
    outliers = final_tumor_vol.loc[(final_tumor_vol < lower_bound) | (final_tumor_vol > upper_bound)]
    
    # Append the current drug's tumor volume data to the list
    tumor_vol_list.append(final_tumor_vol)

    # Print potential outliers for the current drug regimen
    print(f"{drug}'s potential outliers: {outliers}")



# In[20]:


# Generate a box plot that shows the distribution of the tumor volume for each treatment group
orange_out = dict(markerfacecolor="red", markersize=2)
plt.figure(figsize=(8, 6))
plt.boxplot(tumor_vol_list, labels=treatment_list, flierprops=orange_out)
plt.title('Tumor Volume Distribution by Treatment Regimen')
plt.ylabel('Final Tumor Volume (mm3)')
plt.xlabel('Drug Regimen')
plt.grid(True)
plt.show()


# ## Line and Scatter Plots

# In[21]:


# Generate a line plot of tumor volume vs. time point for a single mouse treated with Capomulin
capomulin_line = clean_mouse_study[clean_mouse_study["Drug Regimen"]=="Capomulin"]
mousedata = capomulin_line[capomulin_line["Mouse ID"] == "l509"]
plt.plot(mousedata["Timepoint"], mousedata["Tumor Volume (mm3)"])
plt.xlabel(["Timepoint"])
plt.ylabel("Tumour Volume (mm3)")
plt.title("Capomulin Treatment of mouse l509")
plt.show


# In[22]:


# Step 1: Filter data for Capomulin regimen
capomulin_scatter = clean_mouse_study[clean_mouse_study["Drug Regimen"] == "Capomulin"]

# Step 2: Group by Mouse ID and calculate averages
capomulin_avg = capomulin_scatter.groupby("Mouse ID").mean(numeric_only = True)

# Step 3: Plotting
plt.scatter(capomulin_avg["Weight (g)"], capomulin_avg["Tumor Volume (mm3)"])
plt.xlabel("Weight (g)")
plt.ylabel("Average Tumor Volume (mm3)")
plt.title("Average Tumor Volume vs. Weight for Capomulin Treatment")
plt.show()


# ## Correlation and Regression

# In[28]:


# Calculate the correlation coefficient and a linear regression model
# for mouse weight and average observed tumor volume for the entire Capomulin regimen
corr = st.pearsonr(capomulin_avg["Weight (g)"],capomulin_avg["Tumor Volume (mm3)"])
print(f" The correlation between mouse weight and the average tumor volume is {corr}")


# In[ ]:




