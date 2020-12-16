/**
 * particle_filter.cpp
 *
 * Created on: Dec 12, 2016
 * Author: Tiffany Huang
 */

#include "particle_filter.h"

#include <math.h>
#include <algorithm>
#include <iostream>
#include <iterator>
#include <numeric>
#include <random>
#include <string>
#include <vector>

#include "helper_functions.h"

using std::string;
using std::vector;
using std::normal_distribution;
using std::cout;
using std::endl;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
  /**
   * TODO: Set the number of particles. Initialize all particles to 
   *   first position (based on estimates of x, y, theta and their uncertainties
   *   from GPS) and all weights to 1. 
   * TODO: Add random Gaussian noise to each particle.
   * NOTE: Consult particle_filter.h for more information about this method 
   *   (and others in this file).
   */
  num_particles = 50;  // TODO: Set the number of particles
  Particle my_particle;

  // Create normal distribution for x, y and thtea
  std::default_random_engine gen;
  normal_distribution<double> dist_x(x, std[0]);
  normal_distribution<double> dist_y(y, std[1]);
  normal_distribution<double> dist_theta(theta, std[2]);


  for (int i = 0; i < num_particles; ++i) {
    // Sample from x, y and theta from previous distributions
    my_particle.id = i;
    my_particle.x = dist_x(gen);
    my_particle.y = dist_y(gen);
    my_particle.theta = dist_theta(gen);
    my_particle.weight = 1;


    particles.push_back(my_particle);

    cout << "id = " << my_particle.id << " x = " << my_particle.x << " y = " <<  my_particle.y 
        << " theta = " << my_particle.theta << " w = " << my_particle.weight << endl;
 
  }

  is_initialized = true;

}

void ParticleFilter::prediction(double delta_t, double std_pos[], 
                                double velocity, double yaw_rate) {
  /**
   * TODO: Add measurements to each particle and add random Gaussian noise.
   * NOTE: When adding noise you may find std::normal_distribution 
   *   and std::default_random_engine useful.
   *  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
   *  http://www.cplusplus.com/reference/random/default_random_engine/
   */
  // Variables to store previous coordinates
  double x_0, y_0, theta_0;

  // Create noise process for x, y and thtea
  std::default_random_engine gen;
  normal_distribution<double> noise_x(0, std_pos[0]);
  normal_distribution<double> noise_y(0, std_pos[1]);
  normal_distribution<double> noise_theta(0, std_pos[2]);

  for (int i = 0; i < num_particles; ++i) {

    x_0 = particles[i].x;
    y_0 = particles[i].y;
    theta_0 = particles[i].theta;

    particles[i].x = x_0 + ( sin(theta_0 + yaw_rate*delta_t ) - sin(theta_0) )*velocity/yaw_rate + noise_x(gen);
    particles[i].y = y_0 + ( cos(theta_0) - cos(theta_0 + yaw_rate*delta_t ) )*velocity/yaw_rate + noise_y(gen);
    particles[i].theta = theta_0 + yaw_rate*delta_t + noise_theta(gen);    
    cout << "Prediction - id = " << particles[i].id << velocity << endl;
  }
  

}

void ParticleFilter::dataAssociation(vector<LandmarkObs> predicted, 
                                     vector<LandmarkObs>& observations) {
  /**
   * TODO: Find the predicted measurement that is closest to each 
   *   observed measurement and assign the observed measurement to this 
   *   particular landmark.
   * NOTE: this method will NOT be called by the grading code. But you will 
   *   probably find it useful to implement this method and use it as a helper 
   *   during the updateWeights phase.
   */
  
  // for each measurment/observation  pick the nearest predicted measurment 
  // then assign that landmark to the measurment

  for (uint obs_index = 1; obs_index < observations.size(); obs_index++) {

    double minimum_dist = dist( observations[obs_index].x, observations[obs_index].y, predicted[0].x, predicted[0].y );
    int minimum_dist_index = 0;

    for (uint i = 0; i < predicted.size(); i++) {

      double distance_to_obs = dist( observations[obs_index].x, observations[obs_index].y, predicted[i].x, predicted[i].y );

      if (distance_to_obs < minimum_dist) {
        minimum_dist = distance_to_obs;
        minimum_dist_index = obs_index;
      }

    }

    observations[obs_index].id = predicted[minimum_dist_index].id;

  }  



}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
                                   const vector<LandmarkObs> &observations, 
                                   const Map &map_landmarks) {
  /**
   * TODO: Update the weights of each particle using a mult-variate Gaussian 
   *   distribution. You can read more about this distribution here: 
   *   https://en.wikipedia.org/wiki/Multivariate_normal_distribution
   * NOTE: The observations are given in the VEHICLE'S coordinate system. 
   *   Your particles are located according to the MAP'S coordinate system. 
   *   You will need to transform between the two systems. Keep in mind that
   *   this transformation requires both rotation AND translation (but no scaling).
   *   The following is a good resource for the theory:
   *   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
   *   and the following is a good resource for the actual equation to implement
   *   (look at equation 3.33) http://planning.cs.uiuc.edu/node99.html
   */

  double sigma_xx = std_landmark[0]*std_landmark[0];
  double sigma_yy = std_landmark[1]*std_landmark[1];
  double normalization_factor = 1 / (2 * M_PI * std_landmark[0] * std_landmark[1]);


  int total_landmarks = map_landmarks.landmark_list.size();

  vector<LandmarkObs> observations_local = observations;
  // For each particule find landmarks that are within the sensor_range
  // and predict measurements corresponding to those landmarks
  vector< vector<LandmarkObs> > predicted_measurements(num_particles);

  for (int particle_index = 0; particle_index < num_particles; particle_index++) {

    int landmarks_in_range = 0;
    double distance_to_landmark;

    for ( int landmark_index = 0; landmark_index < total_landmarks; landmark_index++ ) {

      double landmark_x = map_landmarks.landmark_list[landmark_index].x_f;
      double landmark_y = map_landmarks.landmark_list[landmark_index].y_f;

      distance_to_landmark = dist( particles[particle_index].x, particles[particle_index].y, landmark_x, landmark_y );

      if ( distance_to_landmark <= sensor_range ) {
        landmarks_in_range++;
        particles[particle_index].associations.push_back(landmark_index);
        
        // convert landmark coordinates from global frame to vehicle frame
        double landmark_x_local = landmark_x*cos(particles[particle_index].theta) + landmark_y*sin(particles[particle_index].theta) - particles[particle_index].x;
        double landmark_y_local = -landmark_x*sin(particles[particle_index].theta) + landmark_y*cos(particles[particle_index].theta) - particles[particle_index].y;

        LandmarkObs landmark_obs = {landmark_index, landmark_x_local, landmark_y_local };

        predicted_measurements[particle_index].push_back(landmark_obs);
      }


      
    }
    // All landmarks that could be seen by the current particule are added
    // predicted measurments for these landmars were computed
    // moving to data association ( predicted measurments <-> actual observations )
    dataAssociation( predicted_measurements[particle_index], observations_local );

    // here observations_local are assigned to their corresponding landmark
    // update weights for current particle
    // particles[particle_index].weight = 1;
    double weight = 1.0;

    for (uint obs_index = 0; obs_index < observations_local.size();  obs_index++) {

      double mu_x;
      double mu_y;

      // find predicted measurment corresponding to the current observation
      for (uint i = 0; i < predicted_measurements[particle_index].size(); i++) {

        if (observations_local[obs_index].id == predicted_measurements[particle_index][i].id) {
          mu_x = predicted_measurements[particle_index][0].x;
          mu_y = predicted_measurements[particle_index][0].y;
          break;
        }
      }

      double dx_2 = pow( observations_local[obs_index].x - mu_x , 2 );
      double dy_2 = pow( observations_local[obs_index].y - mu_y , 2 );

      weight *= exp( -1/2*( (sigma_xx*dx_2) + (sigma_yy*dy_2) ) ) * normalization_factor;

    }

    particles[particle_index].weight = weight;
    
    // DEBUG
    cout << "Landmarks in range for particule " << particle_index << " => " << landmarks_in_range << endl;


  }

  cout << "Observation vector size => " << observations.size() <<
  " - Map landmarks vector size => " << total_landmarks << endl;
 
  // compute predicted measurement for each particule ( for each landmark)
  // use this pridicted measurment in the dataAssociation phase




  // updating the weight of each particle
  // for each particule we multiply for each measurement 
  for ( int i = 0; i < num_particles; i++ ) {

    particles[i].weight = 1;

  }
  
}

void ParticleFilter::resample() {
  /**
   * TODO: Resample particles with replacement with probability proportional 
   *   to their weight. 
   * NOTE: You may find std::discrete_distribution helpful here.
   *   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
   */


}

void ParticleFilter::SetAssociations(Particle& particle, 
                                     const vector<int>& associations, 
                                     const vector<double>& sense_x, 
                                     const vector<double>& sense_y) {
  // particle: the particle to which assign each listed association, 
  //   and association's (x,y) world coordinates mapping
  // associations: The landmark id that goes along with each listed association
  // sense_x: the associations x mapping already converted to world coordinates
  // sense_y: the associations y mapping already converted to world coordinates
  particle.associations= associations;
  particle.sense_x = sense_x;
  particle.sense_y = sense_y;
}

string ParticleFilter::getAssociations(Particle best) {
  vector<int> v = best.associations;
  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<int>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}

string ParticleFilter::getSenseCoord(Particle best, string coord) {
  vector<double> v;

  if (coord == "X") {
    v = best.sense_x;
  } else {
    v = best.sense_y;
  }

  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<float>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}