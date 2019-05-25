from sklearn.ensemble import RandomForestClassifier
from sklearn.datasets import make_classification
import pandas as pd
import numpy as np

file1= input("File name?"'\n')
#load data set
features1 = pd.read_csv(file1)

#features1
features1=features1.set_index("Unnamed: 0")
#remove row name header
features1.index.name = None

#SAVE COLUMN NAMES
feature1_list = list(features1.columns)

#labels1
labels1 = np.array(features1['HDNA'])
features1= np.array(features1.drop('HDNA', axis = 1))


#dimensions the same
features1.shape
labels1.shape

#split data into training and testing sets; for now I have excluded individual B because they died from hanging,
#and the data seems to have different dispersion

# Using Skicit-learn to split data into training and testing sets
from sklearn.model_selection import train_test_split
# Split the data into training and testing sets
train_features1, test_features1, train_labels1, test_labels1 = train_test_split(features1, labels1, test_size = 0.25, random_state = 42)
#following online tutorial
print('Training Features Shape:', train_features1.shape)
print('Training Labels Shape:', train_labels1.shape)
print('Testing Features Shape:', test_features1.shape)
print('Testing Labels Shape:', test_labels1.shape)

# Import the model we are using
from sklearn.ensemble import RandomForestRegressor
# Instantiate model with 1000 decision trees
rf1 = RandomForestRegressor(n_estimators = 1000, random_state = 42, oob_score= "TRUE")
# Train the model on training data
rf1.fit(train_features1, train_labels1)


# Use the forest's predict method on the test data
predictions = rf1.predict(test_features1)
# Calculate the absolute errors
errors = abs(predictions - test_labels1)


#feature1_list = list(features1.columns)

print('Mean Absolute Error:', round(np.mean(errors), 2))

# Calculate mean absolute percentage error (MAPE)
mape = 100 * (errors / test_labels1)
print('mean absolute percent error:', mape)
#get rid of 0 instance
mape=mape[~np.isinf(mape)]
# Calculate and display accuracy
accuracy = 100 - np.mean(mape)
print('Accuracy:', round(accuracy, 2), '%.')

importances = list(rf1.feature_importances_)
feature1_importances = [(feature, round(importance, 2)) for feature, importance in zip(feature1_list, importances)]
feature1_importances


feature1_importances = sorted(feature1_importances, key = lambda x: x[1], reverse = True)
#[print('Variable: {:20} Importance: {}'.format(*pair)) for pair in feature1_importances];



#combine and plot predictions and test_labels1
df = pd.DataFrame({'x':test_labels1, 'y':predictions})
df.to_csv('test' + file1)

#save important OTUs
important=pd.DataFrame(feature1_importances)
important.to_csv('impOTU' + file1)
#train model using top most important indicators
