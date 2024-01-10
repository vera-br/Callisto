# import packages
from sklearn.datasets import load_iris
from sklearn.model_selection import train_test_split
from xgboost import XGBClassifier

# load data
galileo_wrt_callisto_cphio, B_PDSs = get_pds_data()
galileo_wrt_callisto_cphio_CA = get_closest_approach_data("galileo", "callisto", "cphio", "G")
galileo_wrt_jupiter_SIII = Galileo_trajectories_SIII_from_CPhiO()

# specify orbit
flyby_n = 2

orbit_cphio = galileo_wrt_callisto_cphio["orbit%s" % (flyby_n)]
orbit_SIII = galileo_wrt_jupiter_SIII["orbit%s" % (flyby_n)]
orbit_CA = galileo_wrt_callisto_cphio_CA["CA_orbit%s" % (flyby_n)]
B_PDS = B_PDSs['bfield%s' % (flyby_n)]
B_mag = np.sqrt(B_PDS[1]**2 + B_PDS[2]**2 + B_PDS[3]**2)

time = orbit_cphio[0]

# split into training and test set
X_train, X_test, y_train, y_test = train_test_split(time, B_mag, test_size=0.2, random_state=42)

# Define your hyperparameters
n_estimators_values = [100, 200]
max_depth_values = [None, 5, 10]
learning_rate_values = [0.01, 0.1, 0.2]

# Initialize variables to store the best hyperparameters and performance metrics
best_params = None
best_accuracy = 0 # accuracy is better performance indicator for classification tasks

y_train = y_train.astype(int)
y_test = y_test.astype(int)

# Loop over hyperparameters
for max_depth in max_depth_values:
    for learning_rate in learning_rate_values:
        for n_estimators in n_estimators_values:

            # Define model
            model = XGBClassifier(n_estimators=n_estimators,
                                  max_depth=max_depth,
                                  learning_rate=learning_rate,
                                  objective='binary:logistic',
                                  random_state=42)
            
            # Fit model
            model.fit(X_train, y_train)

            # Calculate prediction
            y_pred = model.predict(X_test)

            # Evaluate performance
            accuracy = accuracy_score(y_test, y_pred)

            # Save parameters if performance is better than previous iterations
            if accuracy > best_accuracy:
                best_accuracy = accuracy
                best_params = {"n_estimators" : n_estimators,
                               "max_depth" : max_depth,
                               "learning_rate" : learning_rate
                               }

# Print best hyperparameters and performance metrics
print("Best Hyperparameters:")
print(best_params)
print("Best accuracy:", best_accuracy)