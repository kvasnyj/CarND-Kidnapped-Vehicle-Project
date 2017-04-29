/*
 * particle_filter.cpp
 *
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
 */

#include <random>
#include <algorithm>
#include <iostream>


#include "particle_filter.h"

std::default_random_engine gen;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
    std::normal_distribution<double> dist_x(x, std[0]);
    std::normal_distribution<double> dist_y(y, std[1]);
    std::normal_distribution<double> dist_psi(theta, std[2]);

    num_particles = 1000;
    for (int i = 0; i < num_particles; i++) {
        Particle p =
                {
                        i, // id
                        dist_x(gen), // x
                        dist_y(gen), // y
                        dist_psi(gen), // theta
                        1 //weight
                };

        particles.push_back(p);
        weights.push_back(1);
    }

    is_initialized = true;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
    for (int i = 0; i < num_particles; i++) {
        double theta =  particles[i].theta;

        double delta_x = 0, delta_y = 0, delta_theta = 0;

        if (yaw_rate == 0) // going straight
        {
            delta_x = velocity * delta_t * sin(theta);
            delta_y = velocity * delta_t * cos(theta);
            delta_theta = 0;
        } else {
            delta_theta = yaw_rate * delta_t;
            delta_x = velocity / yaw_rate * (sin(theta + delta_theta) - sin(theta));
            delta_y = velocity / yaw_rate * (cos(theta) - cos(theta + delta_theta));
        }

        std::normal_distribution<double> dist_x(particles[i].x + delta_x, std_pos[0]);
        std::normal_distribution<double> dist_y(particles[i].y + delta_y, std_pos[1]);
        std::normal_distribution<double> dist_psi(particles[i].theta + delta_theta, std_pos[2]);

        particles[i].x = dist_x(gen);
        particles[i].y = dist_y(gen);
        particles[i].theta = dist_psi(gen);

    }

}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs> &observations) {
    for (int i = 0; i < observations.size(); i++) {
        double x1 = observations[i].x;
        double y1 = observations[i].y;

        double min_d = -1;
        for (int j = 0; j < predicted.size(); j++) {
            double x2 = predicted[j].x;
            double y2 = predicted[j].y;

            double d = dist(x1, y1, x2, y2);
            if (min_d > d || min_d < 0) {
                observations[i].id = j;
                min_d = d;
            }
        }
    }

}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[],
                                   std::vector<LandmarkObs> observations, Map map_landmarks) {
    for (int i = 0; i < num_particles; i++) {
        double x1 = particles[i].x;
        double y1 = particles[i].y;
        double theta = particles[i].theta;

        std::vector<LandmarkObs> predicted_landmarks;
        for (int j = 0; j < map_landmarks.landmark_list.size(); j++) {
            int id = map_landmarks.landmark_list[j].id_i;
            double x2 = map_landmarks.landmark_list[j].x_f;
            double y2 = map_landmarks.landmark_list[j].y_f;

            double delta_x = x2 - x1;
            double delta_y = y2 - y1;

            double d = dist(x1, y1, x2, y2);
            if (d <= sensor_range) {
                // http://planning.cs.uiuc.edu/node99.html
                x2 = delta_x * cos(theta) - delta_y * sin(theta);
                y2 = delta_x * sin(theta) + delta_y * cos(theta);

                LandmarkObs l = {id, x2, y2};
                predicted_landmarks.push_back(l);
            }
        }

        dataAssociation(predicted_landmarks, observations);

        double w = 1;
        for (int j = 0; j < observations.size(); j++) {
            int id = observations[j].id;
            double x1 = observations[j].x;
            double y1 = observations[j].y;

            double x2 = predicted_landmarks[id].x;
            double y2 = predicted_landmarks[id].y;

            double delta_x = x1 - x2;
            double delta_y = y1 - y2;

            double exp1 = exp(-0.5 * (pow(delta_x, 2.0) * std_landmark[0] + pow(delta_y, 2.0) * std_landmark[1]));
            double exp2 = sqrt(2.0 * M_PI * std_landmark[0] * std_landmark[1]);
            w *= exp1 / exp2;
        }

        weights[i] = w;
        particles[i].weight = w;
    }
}

void ParticleFilter::resample() {
    std::vector<Particle> resampled_particles;
    std::discrete_distribution<int> dist_id(weights.begin(), weights.end());

    for (int i = 0; i < particles.size(); i++) {
        int id = dist_id(gen);
        resampled_particles.push_back(particles[id]);
    }

    particles = resampled_particles;
}

void ParticleFilter::write(std::string filename) {
    // You don't need to modify this file.
    std::ofstream dataFile;
    dataFile.open(filename, std::ios::app);
    for (int i = 0; i < num_particles; ++i) {
        dataFile << particles[i].x << " " << particles[i].y << " " << particles[i].theta << "\n";
    }
    dataFile.close();
}
