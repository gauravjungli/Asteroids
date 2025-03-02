import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
from sklearn.model_selection import train_test_split
from sklearn.linear_model import LinearRegression
from sklearn.metrics import r2_score

# Define the system of ODEs
def odes(X, y, phi):
    if len(y) != 2:
        raise ValueError("Input vector y must have length 2.")
    P, psi = y
    # Adding a small value to avoid division by zero
    epsilon = 1e-8
    gn=5
    gt=0
    denominator_1 = -np.sin(phi) + np.cos(2 * psi)
    denominator_2 = 2 * P*  np.sin(phi) * (np.sin(phi) - np.cos(2 * psi))
    if abs(denominator_1)<epsilon:
        print("Denomintor 1 going to zero")
    if abs(denominator_2)<epsilon:
        print("Denomintor 2 going to zero")
 #   if abs(denominator_1) < epsilon or abs(denominator_2) < epsilon:
 #       return [0, 0]  # Return zero derivatives to avoid divide by zero
    dP_dX = 0.1*(gn*np.cos(2 * psi) - gt*np.sin(2*psi)) / denominator_1 if abs(denominator_1)>epsilon else 0

    dpsi_dX = 0.1*( gn*np.sin(phi) * np.sin(2 * psi)-gt*(1-np.sin(phi)*np.cos(2*psi))) / denominator_2 if abs(denominator_2)>epsilon else 0
    
        
    return [dP_dX, dpsi_dX]

# Initial conditions
P0 = 0.1  # Initial value of P
psi0 = np.pi/16  # Initial value of psi
phi = np.pi / 4  # Example value for phi (can be adjusted)

# Time range for the solution
X_span = (0, 10)  # Range of X values
X_eval = np.linspace(X_span[0], X_span[1], 500)  # Points at which to store the solution

# Solve the ODEs
try:
    solution = solve_ivp(odes, X_span, [P0, psi0], args=(phi,), t_eval=X_eval, method='RK45')
    # Extract the solution
    if solution.success:
        X = solution.t
        P = solution.y[0]
        psi = solution.y[1]

        # Plot the results
        plt.figure(figsize=(12, 6))

        plt.subplot(2, 1, 1)
        #plt.plot(X, P, label='P(X)', color='b')
        plt.xlabel('X')
        plt.ylabel('P')
        plt.title('Solution for P(X)')
        plt.grid()
        # Splitting the data into training and testing sets
        X = X.reshape(-1, 1)
        X_train, X_test, y_train, y_test = train_test_split(X, P, test_size=0.2, random_state=42)

        # Creating a linear regression model
        model = LinearRegression()

        # Training the model
        model.fit(X_train, y_train)

        # Predicting using the model
        y_pred = model.predict(X_test)

        # Printing coefficients
        print("Coefficient (Slope):", model.coef_[0])
        print("Intercept:", model.intercept_)
        
        # Calculating and printing R^2 score
        r2 = r2_score(y_test, y_pred)
        print("R^2 Score:", r2)
        # Plotting the results
        plt.plot(X, P, color='blue', label='Original Data')
        plt.plot(X_test, y_pred, color='red', linewidth=2, label='Linear Fit')
        plt.legend()
        #plt.show()

        plt.subplot(2, 1, 2)
        plt.plot(X, (psi), label='psi(X)', color='r')
        plt.xlabel('X',fontsize=16)
        plt.ylabel('psi')
        plt.title('Solution for psi(X)')
        plt.grid()
        plt.legend()
        plt.yticks(fontsize=14)
        plt.xticks(fontsize=14)
        plt.tight_layout()
        plt.show()
    else:
        print("The solver was not successful. Reason:", solution.message)
except ValueError as e:
    print("Error:", e)










