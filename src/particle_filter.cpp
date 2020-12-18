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
using namespace std;


#define EPSILON 0.00001 // Just a small number


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

  // Create normal distribution for x, y and thtea
  std::default_random_engine gen;
  normal_distribution<double> dist_x(x, std[0]);
  normal_distribution<double> dist_y(y, std[1]);
  normal_distribution<double> dist_theta(theta, std[2]);


  for (int i = 0; i < num_particles; ++i) {
    // Sample from x, y and theta from previous distributions

    Particle my_particle;
    my_particle.id = i;
    my_particle.x = x + dist_x(gen);
    my_particle.y = y + dist_y(gen);
    my_particle.theta = theta + dist_theta(gen);
    my_particle.weight = 1.0;


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

  if ( fabs(yaw_rate) < EPSILON ) {
    cout << "Watch out! yaw rate is so small!! - " << yaw_rate << " - " << fabs(yaw_rate) << endl;
  }
  
  // Create noise process for x, y and thtea
  std::default_random_engine gen;
  normal_distribution<double> noise_x(0, std_pos[0]);
  normal_distribution<double> noise_y(0, std_pos[1]);
  normal_distribution<double> noise_theta(0, std_pos[2]);

  for (int i = 0; i < num_particles; ++i) {

    double x_0 = particles[i].x;
    double y_0 = particles[i].y;
    double theta_0 = particles[i].theta;

    if ( fabs(yaw_rate) < EPSILON ) { // Yaw not changing

      particles[i].x += velocity * delta_t * cos( theta_0 ) + noise_x(gen);
      particles[i].y += velocity * delta_t * sin( theta_0 ) + noise_y(gen);

    } else {

      particles[i].x = x_0 + ( sin(theta_0 + yaw_rate*delta_t ) - sin(theta_0) )*velocity/yaw_rate + noise_x(gen);
      particles[i].y = y_0 + ( cos(theta_0) - cos(theta_0 + yaw_rate*delta_t ) )*velocity/yaw_rate + noise_y(gen);
      particles[i].theta = theta_0 + yaw_rate*delta_t + noise_theta(gen);    

    }
    // DEBUG
    //cout << "Prediction - id = " << particles[i].id << "  - x = " << particles[i].x 
    //     << "  - y = " <<  particles[i].y << " - theta = " << particles[i].theta << endl;
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

    LandmarkObs observation = observations[obs_index];

    double minimum_dist = std::numeric_limits<double>::max(); // Intialize with a large number
    int mapId = -1;

    for (uint i = 0; i < predicted.size(); i++) {

      LandmarkObs predicted_lm = predicted[i];

      double distance_to_obs = dist( observation.x, observation.y, predicted_lm.x, predicted_lm.y );

      if (distance_to_obs < minimum_dist) {

        minimum_dist = distance_to_obs;
        mapId = predicted[i].id;

      }

    }

    observations[obs_index].id = mapId;


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


  // This vector will have predicted measurements corresponding to landmarks within the sensor_range
  // for each particule
  vector< vector<LandmarkObs> > predicted_measurements(num_particles);

  for (int particle_index = 0; particle_index < num_particles; particle_index++) {

    

    for ( int landmark_index = 0; landmark_index < total_landmarks; landmark_index++ ) {

      double landmark_x = map_landmarks.landmark_list[landmark_index].x_f;
      double landmark_y = map_landmarks.landmark_list[landmark_index].y_f;
      int landmark_id = map_landmarks.landmark_list[landmark_index].id_i;

      // Calculate distance to landmark
      double distance_to_landmark = dist( particles[particle_index].x, particles[particle_index].y, landmark_x, landmark_y );

      if ( distance_to_landmark <= sensor_range ) {

        LandmarkObs landmark_obs = {landmark_id, landmark_x, landmark_y };

        predicted_measurements[particle_index].push_back(landmark_obs);
      }
      
    }

    // Convert observations to global coordinates
    vector<LandmarkObs> observations_global;

    for (uint obs_index = 0; obs_index < observations.size(); obs_index++) {

      double observations_global_x = observations[obs_index].x*cos(particles[particle_index].theta) - observations[obs_index].y*sin(particles[particle_index].theta) + particles[particle_index].x;
      double observations_global_y = observations[obs_index].x*sin(particles[particle_index].theta) + observations[obs_index].y*cos(particles[particle_index].theta) + particles[particle_index].y;

      LandmarkObs landmark_obs_global = {observations[obs_index].id, observations_global_x, observations_global_y }; // id is not relevent here as it will be populated after in the data association phase

      observations_global.push_back(landmark_obs_global);

    }  

    // All landmarks that could be seen by the current particule are added
    // predicted measurments for these landmars were computed
    // moving to data association ( predicted measurments <-> actual observations ) // using external so far
    dataAssociation( predicted_measurements[particle_index], observations_global );

    // here observations_local are assigned to their corresponding landmark
    // update weights for current particle
    // particles[particle_index].weight = 1;
    double weight = 1.0; // reset weight

    for (uint obs_index = 0; obs_index < observations_global.size();  obs_index++) {

      double prediction_x;
      double prediction_y;

      int obs_id = observations_global[obs_index].id;
      double obs_x = observations_global[obs_index].x;
      double obs_y = observations_global[obs_index].y;

      // find predicted measurment corresponding to the current observation after data association
      for (uint i = 0; i < predicted_measurements[particle_index].size(); i++) {

        if ( obs_id == predicted_measurements[particle_index][i].id) {

          prediction_x = predicted_measurements[particle_index][i].x;
          prediction_y = predicted_measurements[particle_index][i].y;

          break;
        }
      }
      // DEBUG
      //cout << "Particle id = " << particle_index << " - Observation index = " << obs_index
      //<< " - mu_x = " << mu_x << " - mu_y = " << mu_y << endl;

      double dx_2 = pow( obs_x - prediction_x , 2 );
      double dy_2 = pow( obs_y - prediction_y , 2 );

      double weight_i = exp( -( dx_2 / (2 * sigma_xx) + dy_2 / (2 * sigma_yy) ) ) * normalization_factor;

      // avoid loosing all info from weight in case weight_i == 0
      particles[particle_index].weight  = (weight_i == 0.0 ) ? EPSILON : weight * weight_i;

      // DEBUG
      //cout << "Particle id = " << particle_index << " - Intermediate weight = " << weight_i 
      //<< " - dx_2 = " << dx_2 << " - dy_2 = " << dy_2 << endl;


    }

    // cout << "particle  " << particle_index << " Weight is => " << particles[particle_index].weight << endl;

    
    // DEBUG
    //cout << "Landmarks in range for particule " << particle_index << " => " << landmarks_in_range << endl;


  }

  //cout << "Observation vector size => " << observations.size() <<
  //" - Map landmarks vector size => " << total_landmarks << endl;
 
  // compute predicted measurement for each particule ( for each landmark)
  // use this pridicted measurment in the dataAssociation phase





  
  
}

void ParticleFilter::resample() {
  /**
   * TODO: Resample particles with replacement with probability proportional 
   *   to their weight. 
   * NOTE: You may find std::discrete_distribution helpful here.
   *   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
   */

  // Find the maximum weight a particle has
  vector<double> pf_weights;
  double max_weights = numeric_limits<double>::min(); // just to initalize the maximum finding

  for (int i = 0; i < num_particles; i++) {

    pf_weights.push_back(particles[i].weight);
    
    if ( particles[i].weight > max_weights ) {
      max_weights = particles[i].weight;
    }

  }

  std::default_random_engine gen;

  // Create uniform distributions.
  uniform_real_distribution<double> double_dist(0.0, max_weights);
  uniform_int_distribution<int> int_dist(0, num_particles - 1);

  // Generate a random index.
  int index = int_dist(gen);

  double beta = 0.0;

  // Resampling wheel
  vector<Particle> resampledParticles;

  for(int i = 0; i < num_particles; i++) {

    beta += double_dist(gen) * 2.0;

    while( beta > pf_weights[index]) {
      beta -= pf_weights[index];
      index = (index + 1) % num_particles;
    }

    resampledParticles.push_back(particles[index]);

  }

  particles = resampledParticles;

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