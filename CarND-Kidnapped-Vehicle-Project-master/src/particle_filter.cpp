/*
 * particle_filter.cpp
 *
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
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

#include "particle_filter.h"

using namespace std;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).
num_particles=100;
	default_random_engine gen;
	normal_distribution<double> N_x(x, std[0]);
	normal_distribution<double> N_y(y, std[1]);
	normal_distribution<double> N_theta(theta, std[2]);
	for(int i=0; i<num_particles; i++)
	{
		Particle particle;
		particle.id =i;
		particle.x= N_x(gen);
		particle.y= N_y(gen);
		particle.theta= N_theta(gen);
		particle.weight=1;

		particles.push_back(particle);
		weights.push_back(1);
	}
	is_initialized=true;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/
default_random_engine gen;
	for (int i=0;i<num_particles;i++)
	{
		double new_x;
		double new_y;
		double new_theta;
		if (yaw_rate==0)
		{
			new_x=particles[i].x+velocity*delta_t*cos(particles[i].theta);
			new_y=particles[i].y+velocity*delta_t*sin(particles[i].theta);
			new_theta=particles[i].theta;
		}
		else
		{
			new_x=particles[i].x+(velocity/yaw_rate)*(sin(particles[i].theta+yaw_rate*delta_t)- sin(particles[i].theta));
			new_y=particles[i].y+(velocity/yaw_rate)*(cos(particles[i].theta)- cos(particles[i].theta+yaw_rate*delta_t));
			new_theta= particles[i].theta+yaw_rate*delta_t;
		}
		normal_distribution<double> N_x(new_x, std_pos[0]);
		normal_distribution<double> N_y(new_y, std_pos[1]);
		normal_distribution<double> N_theta(new_theta, std_pos[2]);
		particles[i].x=N_x(gen);
		particles[i].y=N_y(gen);
		particles[i].theta=N_theta(gen);
	}
}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.
	for (int i = 0; i < observations.size(); i++) {
		double shortest_distance = 50.0;
		int closest_landmark_id = -1;
		for (int k = 0; k < predicted.size(); k++) {
		  double current_dist = dist(observations[i].x, observations[i].y, predicted[k].x, predicted[k].y);

		  if (current_dist < shortest_distance) {
		    shortest_distance = current_dist;
		    closest_landmark_id = predicted[k].id;
		  }
		}
		observations[i].id = closest_landmark_id;
	}
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
		const std::vector<LandmarkObs> &observations, const Map &map_landmarks) {
	// TODO: Update the weights of each particle using a mult-variate Gaussian distribution. You can read
	//   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
	// NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
	//   according to the MAP'S coordinate system. You will need to transform between the two systems.
	//   Keep in mind that this transformation requires both rotation AND translation (but no scaling).
	//   The following is a good resource for the theory:
	//   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
	//   and the following is a good resource for the actual equation to implement (look at equation 
	//   3.33
	//   http://planning.cs.uiuc.edu/node99.html
  double weight_normalizer = 0.0;

  for (int i = 0; i < num_particles; i++) {
    double particle_x = particles[i].x;
    double particle_y = particles[i].y;
    double particle_theta = particles[i].theta;

    //Transfome in map coordinates
    vector<LandmarkObs> observations_on_map_coor;

    //Transform observations from vehicle's co-ordinates to map co-ordinates.
    for (int k = 0; k < observations.size(); k++) {
      LandmarkObs observations_map;
      observations_map.id = k;
      observations_map.x = particles[i].x + (cos(particles[i].theta) * observations[k].x) - (sin(particles[i].theta) * observations[k].y);
      observations_map.y = particles[i].y + (sin(particles[i].theta) * observations[k].x) + (cos(particles[i].theta) * observations[k].y);
      observations_on_map_coor.push_back(observations_map);
    }

    // only look on the landmarks that are within the sensor range
    vector<LandmarkObs> predicted;
    for (int k = 0; k < map_landmarks.landmark_list.size(); k++) {
      if ((fabs((particle_x - map_landmarks.landmark_list[k].x_f)) <= sensor_range) && (fabs((particle_y - map_landmarks.landmark_list[k].y_f)) <= sensor_range)) {
        predicted.push_back(LandmarkObs {map_landmarks.landmark_list[k].id_i, map_landmarks.landmark_list[k].x_f, map_landmarks.landmark_list[k].y_f});
      }
    }

    // Associate each observation to a landmark
    dataAssociation(predicted, observations_on_map_coor);

    // update the wieghts
    particles[i].weight = 1.0;

    double normalizer = (1.0/(2.0 * M_PI * std_landmark[0] * std_landmark[1]));
    
    /*Calculate the weight of particle based on the multivariate Gaussian probability function*/
    for (int k = 0; k < observations_on_map_coor.size(); k++) {
      particles[i].weight= 1.0;
      for (int m = 0; m < predicted.size(); m++) {
        if (observations_on_map_coor[k].id == predicted[m].id) {
          particles[i].weight*= normalizer *\
           exp(-1.0 * ((pow((observations_on_map_coor[k].x - predicted[m].x), 2)/(2.0 * std_landmark[0]*std_landmark[0])) + \
           	(pow((observations_on_map_coor[k].y - predicted[m].y), 2)/(2.0 * std_landmark[1]*std_landmark[1]))));
        }
      }
    }
    weight_normalizer += particles[i].weight;
  }

  // normalize the wieghts
  for (int i = 0; i < particles.size(); i++) {
    particles[i].weight /= weight_normalizer;
    weights[i] = particles[i].weight;
  }
}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
default_random_engine gen;
	discrete_distribution<int> distribution(weights.begin(),weights.end());
	vector<Particle> resample_particles;
	for(int i=0;i<num_particles;i++)
	{
		resample_particles.push_back(particles[distribution(gen)]);
	}
	particles=resample_particles;
}

Particle ParticleFilter::SetAssociations(Particle& particle, const std::vector<int>& associations, 
                                     const std::vector<double>& sense_x, const std::vector<double>& sense_y)
{
    //particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
    // associations: The landmark id that goes along with each listed association
    // sense_x: the associations x mapping already converted to world coordinates
    // sense_y: the associations y mapping already converted to world coordinates

    particle.associations= associations;
    particle.sense_x = sense_x;
    particle.sense_y = sense_y;
}

string ParticleFilter::getAssociations(Particle best)
{
	vector<int> v = best.associations;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<int>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseX(Particle best)
{
	vector<double> v = best.sense_x;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseY(Particle best)
{
	vector<double> v = best.sense_y;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
