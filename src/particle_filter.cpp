/*
 * particle_filter.cpp
 *
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
 *  
 *  Completed on: May 27, 2017
 *      Author: Sebastián Lucas Sampayo
 */

#include <random>
#include <algorithm>
#include <iostream>
#include <numeric>
#include <math.h> 
#include <iostream>
#include <sstream>
#include <string>
#include <iterator>

#include "helper_functions.h"
#include "map.h"
#include "particle_filter.h"

// #define DEBUG

/**
 *  My previous bugs:
 *    - x error keeps increasing: Solution: 
 *          use sensor_range to take into account only closest landmarks predictions.
 *    - yaw error dramatically increases when making a turn: Solution: 
 *          remove the normalization, the simulator is probably expecting theta in the range [0, 2pi]
 */

using namespace std;

// ------------------------------------------------------------------------------------------------
// -- Utility functions, prototypes
// ------------------------------------------------------------------------------------------------

/**
 *  Normalize an angle between [-pi, +pi]
 */
void normalizeAngle(double& angle);

/**
 *  Transform landmark observations from vehicle coordinates to global coordinates.
 *  @param observations Vector of landmark observations
 *  @param particle Particle object containing the state vector
 */
void transformObservations(std::vector<LandmarkObs>& observations, Particle particle);

/**
 *  Transform landmark observations from vehicle coordinates to global coordinates.
 *  @param[in] x Point X coordinate in vehicle system.
 *  @param particle Particle object containing the state vector
 */
// void transformPoint(double& x_g, double& y_g, const double x_v, const double y_v, Particle particle);

/**
 *  Calculate gaussian pdf in point x
 *  @param x Point to calculate pdf
 *  @param mean Mean of the Gaussian distribution
 *  @param std Standard deviation of the Gaussian distribution
 */
double gaussianPdf(double x, double mean, double std);

/**
 *  Debug: Print all particles
 */
void printParticles(vector<Particle> particles, vector<double> weights);

/**
 *  Debug: Print a vector of landmarks
 */
void printLandmarks(std::vector<LandmarkObs> landmarks, const char* title);
void printLandmarks(Map map);

// ------------------------------------------------------------------------------------------------
// -- ParticleFilter methods
// ------------------------------------------------------------------------------------------------

// ------------------------------------------------------------------------------------------------
void ParticleFilter::init(double x, double y, double theta, double std[]) {
  // Set the number of particles. Initialize all particles to first position (based on estimates of 
  //   x, y, theta and their uncertainties from GPS) and all weights to 1. 
  // Add random Gaussian noise to each particle.
  // NOTE: Consult particle_filter.h for more information about this method (and others in this file).

#ifdef DEBUG
  cout << "init()" << endl;
#endif

  // Less than 100 particles gives worst errors: not enough sampling
  // More than 1000 particles gives worst errors too: probably a performance issue, 
  //                                                  taking too long to run the algorithm.
  // I have found that 100 particles is enough to get the best errors and performance
  num_particles = 100;

  // Clear ParticleFilter weights
  weights.clear();

  // Random distribution initialization
  // @TODO: this can be improved by creating a global normal distribution (0,1) and making a change
  //    of variable each time we need it, taking into account mean and std dev.
  default_random_engine gen;
  
  const double std_x = std[0];
  const double std_y = std[1];
  const double std_theta = std[2];
  normal_distribution<double> dist_x(x, std_x);
  normal_distribution<double> dist_y(y, std_y);
  normal_distribution<double> dist_theta(theta, std_theta);
  
  // Initialization for each particle
  for (int i = 0; i < num_particles; ++i)
  {
    double weight = 1;
    
    weights.push_back(weight);
    
    Particle particle;
    particle.id = i;
    particle.x = dist_x(gen);
    particle.y = dist_y(gen);
    particle.theta = dist_theta(gen);
    particle.weight = weight;
    
    particles.push_back(particle);
  }

  is_initialized = true;

#ifdef DEBUG
  printParticles(particles, weights);
#endif
}

// ------------------------------------------------------------------------------------------------
void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
  // Add measurements to each particle and add random Gaussian noise.

#ifdef DEBUG
  cout << "prediction()" << endl;
#endif

  // Random distribution initialization
  // @TODO: this can be improved by creating a global normal distribution (0,1) and making 
  // a change of variable each time we need it, taking into account mean and std dev.
  default_random_engine gen;
  
  const double std_x = std_pos[0];
  const double std_y = std_pos[1];
  const double std_theta = std_pos[2];
  
  normal_distribution<double> dist_x(0, std_x);
  normal_distribution<double> dist_y(0, std_y);
  normal_distribution<double> dist_theta(0, std_theta);
  
  // Prediction for each particle
  for (Particle& particle : particles )
  {
    double prev_x = particle.x;
    double prev_y = particle.y;
    double prev_theta = particle.theta;
    
    // Motion model
    if (abs(yaw_rate) > 1e-4)
    {
      particle.x += velocity/yaw_rate * (sin(prev_theta + yaw_rate*delta_t) - sin(prev_theta));
      particle.y += velocity/yaw_rate * (cos(prev_theta) - cos(prev_theta + yaw_rate*delta_t));
      particle.theta += yaw_rate * delta_t;
    }
    else
    {
      // Division by zero guard
      particle.x += velocity * delta_t * cos(prev_theta);
      particle.y += velocity * delta_t * sin(prev_theta);
      particle.theta = prev_theta;
    }
    
    // Add noise
    particle.x += dist_x(gen);
    particle.y += dist_y(gen);
    particle.theta += dist_theta(gen);
    
    // Normalization
    // The simulator is probably expecting [0, 2pi] instead of [-pi, pi], so don't normalize
    // normalizeAngle(particle.theta);
  }
  
#ifdef DEBUG
  printParticles(particles, weights);
#endif
}

// ------------------------------------------------------------------------------------------------
void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
  // Find the predicted measurement that is closest to each observed measurement and assign the 
  //   observed measurement to this particular landmark.
  // -> Assigns an id to each observation

#ifdef DEBUG
  cout << "dataAssociation()" << endl;
#endif

  // For each observation, search the closest predicted distance
  for (LandmarkObs& observation : observations)
  {
    double min_dist = 100000; // start with a super high distance
    for (LandmarkObs& pred : predicted)
    {
      double d = dist(pred.x, pred.y, observation.x, observation.y);
      // if it's closest, update minimum distance and id.
      if (d < min_dist)
      {
        min_dist = d;
        observation.id = pred.id;
      }
    }
    // @Optimization: At this point we could remove the closest id from the predicted list because 
    // it's not going to be associated with another observation, it would only add unnecessary operations.
  }
  
#ifdef DEBUG
  // printLandmarks(observations, "Observation");
  // printLandmarks(predicted, "Prediction");
#endif
}

// ------------------------------------------------------------------------------------------------
void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
    std::vector<LandmarkObs> observations, Map map_landmarks) {
  // Update the weights of each particle using a mult-variate Gaussian distribution.
  // NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
  //   according to the MAP'S coordinate system. You will need to transform between the two systems.
  //   Keep in mind that this transformation requires both rotation AND translation (but no scaling).
  //   The following is a good resource for the theory:
  //   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
  //   and the following is a good resource for the actual equation to implement (look at equation 
  //   3.33. Note that you'll need to switch the minus sign in that equation to a plus to account 
  //   for the fact that the map's y-axis actually points downwards.)
  //   http://planning.cs.uiuc.edu/node99.html

#ifdef DEBUG
  cout << "updateWeights()" << endl;
#endif

  const double std_x = std_landmark[0];
  const double std_y = std_landmark[1];

  // Accumulated weights, this is for weight normalization
  double weight_acc = 0;
  
  // Clear ParticleFilter weights
  weights.clear();

  // Update for each particle
  for (Particle& particle : particles)
  {
#ifdef DEBUG
    cout << "-------- Particle id: " << particle.id << 
      " | x: " << particle.x << 
      " | y: " << particle.y << 
      " | theta: " << particle.theta << 
      " | weight: " << particle.weight << 
      endl;
#endif
    // Fill a vector of predicted distances to landmarks for this particle:
    // Note: The predicted landmarks are just the position of the map landmarks we have in our 
    // input map. Although we are only considering landmarks that are closest to the particle 
    // position, in the range defined by sensor_range.
    vector<LandmarkObs> predictions;
    for (auto& landmark : map_landmarks.landmark_list)
    {
      const double lx = landmark.x_f;
      const double ly = landmark.y_f;
      const double distance = dist(lx, ly, particle.x, particle.y);
      
      if (distance <= sensor_range)
      {
        LandmarkObs prediction;
        prediction.id = landmark.id_i;
        prediction.x = lx;
        prediction.y = ly;
        predictions.push_back(prediction);
      }
    }

    // Convert observations from vehicle to global coordinates
    // Create a copy of the observations that will be converted to global coordinates
    std::vector<LandmarkObs> observations_global(observations.size());
    copy(observations.begin(), observations.end(), observations_global.begin());
    transformObservations(observations_global, particle);

    // Associate observations to landmarks, based on predictions for this particle
    dataAssociation(predictions, observations_global);

    // Initialize probability
    double prob = 1;

    // Calculate the posterior probability for each observation
    for (LandmarkObs& observation : observations_global)
    {
      // Search the prediction associated with this landamark id:
      // Optimization: We could make a hash map instead of a lineal search to improve performance
      const int id = observation.id;

      for (LandmarkObs& prediction : predictions)
      {
        if (prediction.id == id)
        {
          // found!
          // For this project we assume that the different axis measurements are independent
          // This allows us to calculate the bivariate distribution as the product of 2 univariate distributions.
          prob *= gaussianPdf(observation.x, prediction.x, std_x);
          prob *= gaussianPdf(observation.y, prediction.y, std_y);
#ifdef DEBUG
          // cout << prob << endl;
#endif
          // done, stop searching
          break;
        }
      }
    }
    
    particle.weight = prob;

    // Accumulate weights
    weight_acc += prob;
  }

  // Normalize weights
  for (Particle& particle : particles)
  {
    particle.weight /= weight_acc;
    weights.push_back(particle.weight);
  }

#ifdef DEBUG
  printParticles(particles, weights);
#endif
}

// ------------------------------------------------------------------------------------------------
void ParticleFilter::resample() {
  // Resample particles with replacement with probability proportional to their weight. 

#ifdef DEBUG
  cout << "resample()" << endl;
#endif
  std::random_device rd;
  std::mt19937 generator_wts(rd());
  std::discrete_distribution<int> d(weights.begin(), weights.end());
  std::vector<Particle> resampledParticles(num_particles);

  // resampledParticles.resize(num_particles);
  for (int i=0; i < num_particles; i++) 
  {
    int idx = d(generator_wts);
    resampledParticles[i] = particles[idx];
  }
  
#ifdef DEBUG
  cout << "resampledParticles.size() = "<< resampledParticles.size() << endl;
#endif

  particles = resampledParticles;

#ifdef DEBUG
  printParticles(particles, weights);
#endif
}

// ------------------------------------------------------------------------------------------------
Particle ParticleFilter::SetAssociations(Particle particle, std::vector<int> associations, std::vector<double> sense_x, std::vector<double> sense_y)
{
  //particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
  // associations: The landmark id that goes along with each listed association
  // sense_x: the associations x mapping already converted to world coordinates
  // sense_y: the associations y mapping already converted to world coordinates

  //Clear the previous associations
  particle.associations.clear();
  particle.sense_x.clear();
  particle.sense_y.clear();

  particle.associations = associations;
   particle.sense_x = sense_x;
   particle.sense_y = sense_y;

   return particle;
}

// ------------------------------------------------------------------------------------------------
string ParticleFilter::getAssociations(Particle best)
{
  vector<int> v = best.associations;
  stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<int>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}

// ------------------------------------------------------------------------------------------------
string ParticleFilter::getSenseX(Particle best)
{
  vector<double> v = best.sense_x;
  stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}

// ------------------------------------------------------------------------------------------------
string ParticleFilter::getSenseY(Particle best)
{
  vector<double> v = best.sense_y;
  stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}

// ------------------------------------------------------------------------------------------------
// -- Utility functions, prototypes
// ------------------------------------------------------------------------------------------------

// ------------------------------------------------------------------------------------------------
void normalizeAngle(double& angle)
{
  //angle normalization
  while (angle >  M_PI) angle -= 2.*M_PI;
  while (angle < -M_PI) angle += 2.*M_PI;
}

// ------------------------------------------------------------------------------------------------
void transformObservations(std::vector<LandmarkObs>& observations, Particle particle)
{
  // The new observation is in vehicle coords, so we have to transform it to global coords:

  // Transform the landmark position from vehicle to global coordinates.
  // Rotation coord from vehicle to global(y-axis upwards):
    // [ cos(theta) -sin(theta)
    //   sin(theta) cos(theta)]
  // y-axis inversion:
    // [1 0
    //  0 -1]
  // Rotation coord from vehicle to global(y-axis downwards) (i.e. [T]^{GV}(theta) ):
  // T = y_inversion * rot
    // [  cos(theta) -sin(theta)
    //   -sin(theta) -cos(theta)]

  // Considering the axis.png, no need to invert axis.
  // Translation:
  // x_vehicle = x_global + x_particle
  // [x]^G = [T]^{GV}(theta) * [x]^V + [p]^G, where x is the measurement and p the particle position

  const double yaw = particle.theta;
  const double px = particle.x;
  const double py = particle.y;
  
  for (LandmarkObs& observation : observations)
  {
    const double obs_x = observation.x; // vehicle coords
    const double obs_y = observation.y; // vehicle coords
    
    const double transformed_obs_x = cos(yaw) * obs_x - sin(yaw) * obs_y + px;
    const double transformed_obs_y = sin(yaw) * obs_x + cos(yaw) * obs_y + py;
    
    observation.x = transformed_obs_x; // global coords
    observation.y = transformed_obs_y; // global coords
  }
}

// ------------------------------------------------------------------------------------------------
double gaussianPdf(double x, double mean, double std)
{
  double normal_x_2 = (mean - x) * (mean - x);
  double std2 = std*std;
  double gaussian = exp(-normal_x_2 / std2 / 2.0) / sqrt(2.0 * M_PI * std2);
  return gaussian;
}

// ------------------------------------------------------------------------------------------------
void printParticles(vector<Particle> particles, vector<double> weights)
{
  for (Particle particle : particles)
  {
    cout << "Particle id: " << particle.id << 
      " | x: " << particle.x << 
      " | y: " << particle.y << 
      " | theta: " << particle.theta << 
      " | p.weight: " << particle.weight << 
      endl;
  }
  
  cout << "--- weights ---" << endl;
  
  for (double weight : weights)
  {
    cout << "weights[i]: " << weight << endl;
  }
}

// ------------------------------------------------------------------------------------------------
void printLandmarks(std::vector<LandmarkObs> landmarks, const char* title)
{
  for (LandmarkObs landmark : landmarks)
  {
    cout << title << " id: " << landmark.id << 
      " | x: " << landmark.x << 
      " | y: " << landmark.y << 
      endl;
  }
}

// ------------------------------------------------------------------------------------------------
void printLandmarks(Map map)
{
  for (auto& landmark : map.landmark_list)
  {
    cout << "Map Landmark id: " << landmark.id_i << 
      " | x: " << landmark.x_f << 
      " | y: " << landmark.y_f << 
      endl;
  }
}