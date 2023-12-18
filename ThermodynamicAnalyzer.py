import numpy as np
import matplotlib.pyplot as plt

class ThermodynamicsAnalyzer:
    R = 8.314  # Universal gas constant (J/(mol K))

    def __init__(self, initial_temperature, initial_pressure, initial_volume, mass, length, moles=1):
        self.temperature = initial_temperature
        self.pressure = initial_pressure
        self.volume = initial_volume
        self.mass = mass
        self.length = length
        self.moles = moles

        # Lists to store data for plotting
        self.time_points = []
        self.length_points = []
        self.pressure_points = []
        self.volume_points = []

    def calculate_internal_energy(self):
        # Internal Energy (U) = n * Cv * T
        return self.moles * self.get_specific_heat() * self.temperature

    def calculate_entropy(self):
        # Entropy (S) = n * Cv * ln(T) + n * R * ln(V)
        return (
            self.moles * self.get_specific_heat() * np.log(self.temperature) +
            self.moles * self.R * np.log(self.volume)
        )

    def calculate_specific_heat(self):
        # Specific Heat (Cv) = f/2 * R, assuming a diatomic gas with f degrees of freedom
        degrees_of_freedom = 5  # for a diatomic gas
        return degrees_of_freedom / 2 * self.R

    def update_state_euler(self, time_step):
        # Euler method to approximate changes in length, pressure, and volume over time
        delta_temperature = time_step
        delta_length = self.length * self.alpha * delta_temperature  # Euler method for length change
        delta_pressure = self.moles * self.R * delta_temperature / self.volume
        delta_volume = self.volume * self.beta * delta_temperature  # Euler method for volume change

        self.temperature += delta_temperature
        self.length += delta_length
        self.pressure += delta_pressure
        self.volume += delta_volume

        # Store data for plotting
        self.time_points.append(self.temperature)
        self.length_points.append(self.length)
        self.pressure_points.append(self.pressure)
        self.volume_points.append(self.volume)

    def perform_analysis_euler(self, total_time, time_step):
        for _ in range(int(total_time / time_step)):
            self.update_state_euler(time_step)

        internal_energy = self.calculate_internal_energy()
        entropy = self.calculate_entropy()
        specific_heat = self.calculate_specific_heat()

        # Additional thermodynamic calculations
        density = self.mass / self.volume
        pressure_from_ideal_gas_law = self.moles * self.R * self.temperature / self.volume

        analysis_results = {
            "Internal Energy": internal_energy,
            "Entropy": entropy,
            "Specific Heat": specific_heat,
            "Statistics": {
                "Mean": np.mean(self.volume_points),
                "Standard Deviation": np.std(self.volume_points)
            },
            "Additional Thermodynamic Calculations": {
                "Density": density,
                "Pressure (Ideal Gas Law)": pressure_from_ideal_gas_law
            }
        }

        return analysis_results

    def plot_results(self):
        # Plot length, pressure, and volume over time
        plt.figure(figsize=(12, 8))

        plt.subplot(3, 1, 1)
        plt.plot(self.time_points, self.length_points, label='Length')
        plt.xlabel('Time')
        plt.ylabel('Length')
        plt.legend()

        plt.subplot(3, 1, 2)
        plt.plot(self.time_points, self.pressure_points, label='Pressure')
        plt.xlabel('Time')
        plt.ylabel('Pressure')
        plt.legend()

        plt.subplot(3, 1, 3)
        plt.plot(self.time_points, self.volume_points, label='Volume')
        plt.xlabel('Time')
        plt.ylabel('Volume')
        plt.legend()

        plt.tight_layout()
        plt.show()


# Example usage
if __name__ == "__main__":
    # Input parameters (replace with actual values)
    initial_temperature = 300  # in Kelvin
    initial_pressure = 1  # in atm
    initial_volume = 10  # in cubic meters
    mass = 0.1  # in kg
    length = 2  # in meters
    alpha = 0.00001  # coefficient of linear expansion
    beta = 0.00002  # coefficient of volume expansion

    # Create an instance of ThermodynamicsAnalyzer
    analyzer = ThermodynamicsAnalyzer(initial_temperature, initial_pressure, initial_volume, mass, length)

    # Set coefficients of expansion
    analyzer.alpha = alpha
    analyzer.beta = beta

    # Perform thermodynamic analysis using Euler method and plot results
    results = analyzer.perform_analysis_euler(total
