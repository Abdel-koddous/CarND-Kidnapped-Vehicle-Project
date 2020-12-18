## **Kidnapped Vehicle Project**
#### Abdelkoddous Khamsi

[gif1]: ./writeup_resources/particle_filter_in_action.gif "Live Particle Filter"


The goal of this project is to implement a 2d particle filter to localize a vehicle in a map using some landmarks measurements.


The main steps I went through during the implementation of the Particle filter in the `particle_filter.cpp` file are:
* Particle Filter initialization.
* Particle filter prediction step.
* Particles weights update step.
* Particles resampling.


### **1. Particle Filter initialization:**

* I Created `50` particles whose `x`, `y` and `theta` values are drawn from normal distributions with means and variances given by the GPS module of the vehicle. 
* At the intialization step the `weight` variables of all the particles are set to **1**.  

### **2. Particle Filter prediction step:**

I predict the new values of `x`, `y` and `theta` for each particle using the _Bicycle model_ yaw and velocity equations from the _Motion models_ lesson in the `.prediction(double delta_t, double std_pos[], double velocity, double yaw_rate)` method.

### **3. Particles weights update step:**

The was by far the part were I spent most of my time in the project.
* First, all landmarks that could be seen (within sensor range) by the current particule are added, predicted measurments for these landmars were computed.
* Then, I had to perform **Data Association** between every actual observation and the nearest predicted observation for it in order to associate the measurement with a landmark on the map.
* I update the weight of each particule to account for how likely are the actual current observations from the point of view of the every particle. 


### **4. Particles Resampling.:**
In this part I implemented the resampling wheel following the algorithm discussed in the lesson.

### **5. Particle Filter in action**

Here is an overview of the final implementation from the term 2 simulator.

![alt text][gif1]

The full _mp4_ video can be found [here](./writeup_resources/particle_filter_final_run.mp4).