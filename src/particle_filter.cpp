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
#define NUM_PARTICLES 100

using std::string;
using std::vector;

using namespace std;

default_random_engine gen;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
  /**
   * TODO: Set the number of particles. Initialize all particles to 
   *   first position (based on estimates of x, y, theta and their uncertainties
   *   from GPS) and all weights to 1. 
   * TODO: Add random Gaussian noise to each particle.
   * NOTE: Consult particle_filter.h for more information about this method 
   *   (and others in this file).
   */
  num_particles = NUM_PARTICLES;  // TODO: Set the number of particles
  
    normal_distribution<double> dist_x(x, std[0]);
    normal_distribution<double> dist_y(y, std[1]);
    normal_distribution<double> dist_theta(theta, std[2]);
    
    for (int i=0; i < num_particles; i++)
    {
        Particle sample;
        sample.x = dist_x(gen);
        sample.y = dist_y(gen);
        sample.theta = dist_theta(gen);
        sample.weight = 1.0;
        
        particles.push_back(sample);
        ParticleFilter::weights.push_back(1.0);
    }
    
    ParticleFilter::is_initialized = true;

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
      for (int i=0; i < num_particles; i++)
    {
        normal_distribution<double> dist_x(particles[i].x, std_pos[0]);
        normal_distribution<double> dist_y(particles[i].y, std_pos[1]);
        normal_distribution<double> dist_theta(particles[i].theta, std_pos[2]);
        
        double theta_gauss = dist_theta(gen);
        double x_gauss = dist_x(gen);
        double y_gauss = dist_y(gen);
        
        if (yaw_rate != 0)
        {
            particles[i].x = x_gauss + ((velocity/yaw_rate)*
            (sin(theta_gauss + (yaw_rate * delta_t)) - sin(theta_gauss)));

            particles[i].y = y_gauss + ((velocity/yaw_rate)*
            (cos(theta_gauss) - cos(theta_gauss+(yaw_rate*delta_t))));
        }
        else
        {
            particles[i].x = x_gauss + (velocity * cos(theta_gauss) * delta_t);
            particles[i].y = y_gauss + (velocity * sin(theta_gauss) * delta_t);
            
        }
        particles[i].theta = theta_gauss + yaw_rate * delta_t;
        
        dist_x.reset();
        dist_y.reset();
        dist_theta.reset();
        
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

   double delta_x, delta_y, diff;
   double gauss_norm = 1.0 / (2 * M_PI * std_landmark[0] * std_landmark[1]);


   // Iterate through all particles 
    for (int i=0; i < num_particles; i++)
    {
        particles[i].NN_associations = observations;
        vector<LandmarkObs> temp;

        // For each association, transform from vehicle coordinates to map coordinates

        for (unsigned int j=0; j < particles[i].NN_associations.size(); j++)
        {
            double range = sqrt(pow(observations[j].x,2) + pow(observations[j].y,2));
            if (range < sensor_range)
            {
                particles[i].NN_associations[j].x = (observations[j].x * cos(particles[i].theta)) - (observations[j].y * sin(particles[i].theta)) + particles[i].x;
                
                particles[i].NN_associations[j].y = (observations[j].x * sin(particles[i].theta)) + (observations[j].y * cos(particles[i].theta)) + particles[i].y;
              
                temp.push_back(particles[i].NN_associations[j]);
            }
        }
        particles[i].NN_associations = temp;
        
        // Setting Nearest Neighbour for each prediction of the Particle


        for (unsigned int j=0; j < particles[i].NN_associations.size(); j++)
        {
            double min_diff = 10000;
            for (unsigned int k=0; k < map_landmarks.landmark_list.size(); k++)
            {
                delta_x = particles[i].NN_associations[j].x - map_landmarks.landmark_list[k].x_f;
                delta_y = particles[i].NN_associations[j].y - map_landmarks.landmark_list[k].y_f;
                diff = sqrt(pow(delta_x,2) + pow(delta_y,2));
                
                if (diff < min_diff)
                {
                    min_diff = diff;
                    particles[i].NN_associations[j].id = map_landmarks.landmark_list[k].id_i;
                }
            }
        }
        particles[i].weight = 1.0;

        // Calculating and updating weight for each Nearest Neighbor

        for (unsigned int j=0; j < particles[i].NN_associations.size(); j++)
        {
            Map::single_landmark_s landmark =  map_landmarks.landmark_list[particles[i].NN_associations[j].id - 1];
            double diff_x = particles[i].NN_associations[j].x -landmark.x_f;
            double X = pow(diff_x,2) / (2 * pow(std_landmark[0],2));
        
            double diff_y = particles[i].NN_associations[j].y -landmark.y_f;
            double Y = pow(diff_y,2) / (2 * pow(std_landmark[1],2));
        
            double weight_temp = exp(-(X+Y)) * gauss_norm;
            particles[i].weight *= weight_temp;
        }

        ParticleFilter::weights[i] = particles[i].weight;
    }

}

void ParticleFilter::resample() {
  /**
   * TODO: Resample particles with replacement with probability proportional 
   *   to their weight. 
   * NOTE: You may find std::discrete_distribution helpful here.
   *   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
   */
    discrete_distribution<> d(ParticleFilter::weights.begin(), ParticleFilter::weights.end());
    vector<Particle> resampled_particles;
    for (int i=0; i < ParticleFilter::num_particles; i++)
    {
        resampled_particles.push_back(particles[d(gen)]);
    }
    particles = resampled_particles;
    d.reset();

}

void ParticleFilter::SetAssociations(Particle& best) {
  // particle: the particle to which assign each listed association, 
  //   and association's (x,y) world coordinates mapping
  // associations: The landmark id that goes along with each listed association
  // sense_x: the associations x mapping already converted to world coordinates
  // sense_y: the associations y mapping already converted to world coordinates
    vector<int> temp_associations;
    vector<double> temp_sense_x;
    vector<double> temp_sense_y;
    for(unsigned int j=0; j < best.NN_associations.size(); j++)
    {
        temp_associations.push_back(best.NN_associations[j].id);
        temp_sense_x.push_back(best.NN_associations[j].x);
        temp_sense_y.push_back(best.NN_associations[j].y);
    }
    best.associations = temp_associations;
    best.sense_x = temp_sense_x;
    best.sense_y = temp_sense_y;
    
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