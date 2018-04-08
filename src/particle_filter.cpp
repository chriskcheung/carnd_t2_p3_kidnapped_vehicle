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
	
	// set number of particles
	num_particles = 8;

	// standard deviations for x, y, and theta
	//std[0] for x in meter
	//std[1] for y in meter
	//std[2] for yaw in radian
	
	// create normal gaussian distribution for x, y, theta
	normal_distribution<double> dist_x(x, std[0]);
	normal_distribution<double> dist_y(y, std[1]);
	normal_distribution<double> dist_theta(theta, std[2]);
	
	// generate particles by sampling them base on the first gps positions provided with normal gaussian 
	// distribution within their standard deviations
	for (int i=0; i<num_particles; i++){
		default_random_engine gen;
		double sample_x, sample_y, sample_theta;
		
		sample_x = dist_x(gen);
		sample_y = dist_y(gen);
		sample_theta = dist_theta(gen);
		Particle sample = {i, sample_x, sample_y, sample_theta, 1.0};
		particles.push_back(sample);
		weights.push_back(1);
		//cout << "init: x|y|theta=" << sample_x << "|" << sample_y << "|" << sample_theta << endl; 
		//cout << "part: x|y|theta=" << particles[i].x << "|" << particles[i].y << "|" << particles[i].theta << "|" << particles[i].id << endl; 
	}
	

	// set is_initialized to true
	is_initialized = true;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/

	normal_distribution<double> dist_x(0, std_pos[0]);
	normal_distribution<double> dist_y(0, std_pos[1]);
	normal_distribution<double> dist_theta(0, std_pos[2]);
	
	for (int i=0; i<num_particles; i++){
		double x, y, theta;
		x = particles[i].x;
		y = particles[i].y;
		theta = particles[i].theta;
		default_random_engine gen;
		
		if (yaw_rate != 0.0f){
			x = x + velocity/yaw_rate*( sin(theta + yaw_rate*delta_t) - sin(theta));
			y = y + velocity/yaw_rate*(-cos(theta + yaw_rate*delta_t) + cos(theta));
			theta = theta + yaw_rate*delta_t;
		}
		else {
			x = x + velocity*cos(theta)*delta_t;
			y = y + velocity*sin(theta)*delta_t;
			theta = theta;
		}
		particles[i].x = x +dist_x(gen);
		particles[i].y = y+dist_y(gen);
		particles[i].theta = theta+dist_theta(gen);
		//cout << "predict: x|y|theta=" << particles[i].x << "|" << particles[i].y << "|" << particles[i].theta << "|" << particles[i].id << endl; 
	}
}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.
	// find the shortest distance to landmark
	// predicted    - map landmarks
	// observations - transformed observation of landmarks
	for (int j=0; j<observations.size() && predicted.size()>0; j++){
		double shortest;
		shortest = std::numeric_limits<double>::infinity();
		int    shortest_idx;
		for (int k=0; k<predicted.size(); k++){
			double diff = dist(observations[j].x, observations[j].y, predicted[k].x, predicted[k].y);
			if (diff < shortest){
				shortest = diff;
				shortest_idx = k;
			}
			//cout << "j|k=" << j << "|" << k << " shortest, shortest_idx, observations[j].x, observations[j].y, predicted[k].x, predicted[k].y=" << shortest << "," << shortest_idx << "," << observations[j].x <<","<< observations[j].y <<","<< predicted[k].x <<","<< predicted[k].y << endl;
		}
		if (shortest < std::numeric_limits<double>::infinity()){			
			observations[j]=predicted[shortest_idx];
			//cout << "nearest obs=" << observations[j].x << "," << observations[j].y << endl;
		}
		else {
			cout << "no shortest found. Error" << endl;
		}
		// if shortest is less than certain distance, move it from the predicted(landmark list) to reduce computation time
		//if (shortest <= 2.0f) {
		//	predicted.erase(predicted.begin() + shortest_idx);
		//}
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
	
	double stdx2  = pow(std_landmark[0], 2);
	double stdy2  = pow(std_landmark[1], 2);
			
	for (int i=0; i<num_particles; i++){
		// use sensor_range to minimize landmarks that are too far away from particles
		// only consider those that are within the range the sensor can cover
		vector<LandmarkObs> reduced_map;
		
		// TODO: COMMENTED OUT OPTIMIZATION
		/*
		for (int j=0; j<map_landmarks.landmark_list.size(); j++){
			if (dist(particles[i].x, particles[i].y, map_landmarks.landmark_list[j].x_f, map_landmarks.landmark_list[j].y_f) <= sensor_range){
				LandmarkObs lmObs;
				lmObs.id = map_landmarks.landmark_list[j].id_i;
				lmObs.x  = map_landmarks.landmark_list[j].x_f;
				lmObs.y  = map_landmarks.landmark_list[j].y_f;
				reduced_map.push_back(lmObs);
			}
		}
		//cout << "original map_landmarks number=" << map_landmarks.landmark_list.size() << endl;
		//cout << "reduced  map_landmarks number=" << reduced_map.size() << endl;
		//cout << "P[" << i << "]" << endl;
		//for (int i=0; i<reduced_map.size(); i++)
		//	cout << "reduced_map[" << i<< "]x,y=" << reduced_map[i].x << "," << reduced_map[i].y << endl;
		*/
		
		
		
		// transform particle coordinates into map coordinates
		// for every observation from car's coordinates, assume particles shares the same observations, transform every observation to gps map coordinate
		std::vector<LandmarkObs> transformed_obs;
		for (int j=0; j<observations.size(); j++){
			LandmarkObs tmp;
			tmp.x = particles[i].x + observations[j].x*cos(particles[i].theta) - observations[j].y*sin(particles[i].theta);
			tmp.y = particles[i].y + observations[j].y*cos(particles[i].theta) + observations[j].x*sin(particles[i].theta);
			transformed_obs.push_back(tmp);

			//cout << "original  observation     coordinate x,y=" << observations[j].x << "," << observations[j].y << endl;
			//cout << "transform transformed_obs coordinate x,y=" << transformed_obs[j].x << "," << transformed_obs[j].y << endl;
			
			double shortest = std::numeric_limits<double>::infinity();			
			for (int k = 0; k < map_landmarks.landmark_list.size(); k++) {
				double diff = dist(tmp.x, tmp.y, map_landmarks.landmark_list[k].x_f, map_landmarks.landmark_list[k].y_f);
				if (diff < shortest) {
					shortest = diff;
					transformed_obs[j].id = map_landmarks.landmark_list[k].id_i;
				}
			}
		}
		
		
		// associate newly transformed observations to nearest landmark
		// find the shortest distance from transformed_obs to landmark
		//vector<LandmarkObs> associated_obs = transformed_obs;
		//dataAssociation(reduced_map, associated_obs);
		/*
		for (int j=0; j<observations.size(); j++){
			double shortest;
			shortest = std::numeric_limits<double>::infinity();
			int    shortest_idx;
			
			for (int k = 0; i < map_landmarks.landmark_list.size(); k++) {
				double diff = dist(transformed_obs[j].x, transformed_obs[j].y, map_landmarks.landmark_list[k].x, map_landmarks.landmark_list[k].y);
				if (diff < shortest) {
					shortest = diff;
					shortest_idx = 
				}
			}
			
			for (int k=0; k<predicted.size(); k++){
				double diff = dist(observations[j].x, observations[j].y, predicted[k].x, predicted[k].y);
			if (diff < shortest){
				shortest = diff;
				shortest_idx = k;
			}
			//cout << "j|k=" << j << "|" << k << " shortest, shortest_idx, observations[j].x, observations[j].y, predicted[k].x, predicted[k].y=" << shortest << "," << shortest_idx << "," << observations[j].x <<","<< observations[j].y <<","<< predicted[k].x <<","<< predicted[k].y << endl;
		}
		if (shortest < std::numeric_limits<double>::infinity()){			
			observations[j]=predicted[shortest_idx];
			//cout << "nearest obs=" << observations[j].x << "," << observations[j].y << endl;
		}
		else {
			cout << "no shortest found. Error" << endl;
		}
		*/
		
		
		// update weights of all particles
		// using Multivariate-Gaussian's standard deviation algorithm to find the sum of product between particles and associated as weight update
		// x and y are the observations in map coordinates (transformed_obs) and μx and μy are the coordinates of the nearest landmarks (associated_obs)
		double new_weight = 1.0f;
		//cout << "transformed_obs" << transformed_obs.size() << endl;
		for (int j=0; j<transformed_obs.size(); j++){
			//double x = transformed_obs[j].x;
			double x = map_landmarks.landmark_list[transformed_obs[j].id - 1].x_f;
			//double y = transformed_obs[j].y;
			double y = map_landmarks.landmark_list[transformed_obs[j].id - 1].y_f;
			//double ux = associated_obs[j].x;
			//double uy = associated_obs[j].y;
			double ux = transformed_obs[j].x;
			double uy = transformed_obs[j].y;
			//double diffx2 = pow((x-ux),2.0f);
			double diffx2 = (x-ux) * (x - ux);
			//double diffy2 = pow((y-uy),2.0f);
			double diffy2 = (y-uy) * (y - uy);
			double sqrt2pr= (2.0f*M_PI*std_landmark[0]*std_landmark[1]);
			//cout << "x|y|ux|uy|diffx2|diffy2|sqrt2pr=" << x << "|" << y << "|" << ux << "|" << uy << "|" << diffx2 << "|" << diffy2 << "|" << sqrt2pr << endl;
			cout << "exp=" << exp(-(diffx2/(2.0f*stdx2)) - (diffy2/(2.0f*stdy2)))/sqrt2pr << endl;
			new_weight *= exp(-(diffx2/(2.0f*stdx2)) - (diffy2/(2.0f*stdy2)))/sqrt2pr;
		}
		weights[i] = new_weight*1.0f;
		particles[i].weight = new_weight;
		//cout << "new_weight=" << new_weight << endl;
		cout << "sum product of weights[" << i << "]=" << weights[i] << endl;
	}
	
	/*
	// normalized all weights by dividing each weight element with sum of all weights
	double norm_weight = std::accumulate(weights.begin(), weights.end(), 0.0f);
	//cout << "norm_weight=" << norm_weight << endl;
	if (norm_weight > 0) {
		for (int i =0; i<num_particles; i++){
			weights[i] /= norm_weight*1.0f;
			particles[i].weight /= norm_weight;
			//cout << "normalized weights[" << i << "]=" << weights[i] << endl;
		}
	}
	*/
}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
	// using sabastian's resampling wheel concept to resample
	
	//TODO , need to use particles[i].weight instead of weight
/*
	vector<Particle> new_particles;
	auto max_weight = *max_element(weights.begin(), weights.end());
	// randomized to an index to start with
	int    index = rand()%num_particles;
	double beta = 0.0f;
	for (int i=0, j=0; i<num_particles; i++){
		//default_random_engine gen;
		//cout << "gen=" << gen << endl;

		beta += (rand()%3)*max_weight;
		//cout << "beta=" << beta << "  max_weight=" << max_weight << "starting index=" << index << endl;
		while (weights[index] < beta){
			beta -= weights[index];
			index = (index + 1)%num_particles;
			//cout << "while loop " << j++ << " beta=" << beta << "index=" << index << endl; 
		}
		new_particles.push_back(particles[index]);
		//cout << "new_particles=" << new_particles[i].id << "," << new_particles[i].x << "," << new_particles[i].y << "," << new_particles[i].weight << endl;
	}
	particles = new_particles;
*/

	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution

	default_random_engine gen;
	discrete_distribution<> ddist(weights.begin(), weights.end());
	
	std::vector<Particle> new_particles;
	for (int i = 0; i < num_particles; i++) {
		new_particles.push_back(particles[ddist(gen)]);
	}
	particles = new_particles;
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
