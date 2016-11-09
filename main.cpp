#include <cmath>
#include <random>
#include <iostream>
#include <omp.h>

#define POINTS 15000
#define DIM 5
#define CLUSTERS 50
#define MAXITER 1000

using namespace std;
using namespace std::chrono;

class Point
{
private:
	double* coordinates;
public:
	Point(double* coordinates)
	{
		this->coordinates = coordinates;
	}
	virtual double* get_coordinates()
	{
		return this->coordinates;
	}
};

class Cluster : Point
{
private:
	double* center_buffer;
	int points;
public:
	Cluster(double *coordinates) : Point(coordinates)
	{
		this->center_buffer = new double[DIM];
	}
	double* get_center_buffer()
	{
		return this->center_buffer;
	}
	double* get_coordinates()
	{
		return Point::get_coordinates();
	}
	void add_to_buffer(double* values, int size)
	{
		for (int i = 0; i < size; i++)
		{
			this->center_buffer[i] += values[i];
		}
		this->points++;
	}
	void reset_points()
	{
		this->points = 0;
	}
	int get_points()
	{
		return this->points;
	}
};

double distance(double* point1, double* point2, int size)
{
	double sum = 0;
	for (auto i = 0; i < size; i++)
	{
		sum += pow(point1[i] - point2[i], 2);
	}

	//cout << "Distance is : " << sqrt(sum) << endl;

	return sqrt(sum);
}

void single_threaded(Cluster* clusters[], Point* points[])
{
    // Will be set to 1 if any point-cluster pair changes
    int differ = 0;

    // Tells us how many times we rearranged the point-cluster pairs
    int iterations = 0;

    Cluster* closest_centers[POINTS];

    do {
        // Assume no change until now
        differ = 0;

        iterations++;

        // Reset clusters center buffers and nr of points
        for (auto i = 0; i < CLUSTERS; i++)
        {
            fill(clusters[i]->get_center_buffer(), clusters[i]->get_center_buffer() + DIM, 0);
            clusters[i]->reset_points();
        }

        // Get nearest cluster for each point
        for (auto i = 0; i < POINTS; i++)
        {
            double min_distance = 10000;
            Cluster *min_cluster = nullptr;

            for (auto j = 0; j < CLUSTERS; j++)
            {
                double current_distance = distance(
                        points[i]->get_coordinates(),
                        clusters[j]->get_coordinates(),
                        DIM
                );

                if (min_distance > current_distance)
                {
                    min_distance = current_distance;
                    min_cluster = clusters[j];
                }
            }

            if (closest_centers[i] != min_cluster)
            {
                // We changed a cluster-point pair
                differ = 1;
                closest_centers[i] = min_cluster;
            }

            // Update buffer that will be used to recalculate center
            closest_centers[i]->add_to_buffer(points[i]->get_coordinates(), DIM);
        }

        if (!differ)
            break;

        // Get center for each cluster
        for (auto i = 0; i < CLUSTERS; i++)
        {
            if (clusters[i]->get_points() == 0)
            {
                continue;
            }
            for (auto j = 0; j < DIM; j++)
            {
                clusters[i]->get_coordinates()[j] = clusters[i]->get_center_buffer()[j]
                                                    / clusters[i]->get_points();
            }
        }
    } while (differ && iterations < MAXITER);

    cout << endl << "Result: " << endl;
    cout << clusters[0]->get_coordinates()[0] << "  ";
    cout << clusters[0]->get_coordinates()[1];
}

void multi_threaded(Cluster* clusters[], Point* points[])
{
    // Will be set to 1 if any point-cluster pair changes
    int differ = 0;

    // Tells us how many times we rearranged the point-cluster pairs
    int iterations = 0;

    // Save closest center for each point
    Cluster* closest_centers[POINTS];

    // Lock for updating cluster buffers
    omp_lock_t lock_cluster_buffer;
    omp_init_lock(&lock_cluster_buffer);

    do {
        // Assume no change until now
        differ = 0;

        iterations++;

        // Reset clusters center buffers and nr of points
        #pragma omp parallel num_threads(5)
        {
            #pragma omp for
            for (int i = 0; i < CLUSTERS; i++) {
                fill(clusters[i]->get_center_buffer(), clusters[i]->get_center_buffer() + DIM, 0);
                clusters[i]->reset_points();
            }
        }

        // Get nearest cluster for each point
        #pragma omp parallel num_threads(5)
        {
            #pragma omp for
            for (int i = 0; i < POINTS; i++)
            {
                double min_distance = 10000;
                Cluster *min_cluster = nullptr;

                for (int j = 0; j < CLUSTERS; j++)
                {
                    double current_distance = distance(
                            points[i]->get_coordinates(),
                            clusters[j]->get_coordinates(),
                            DIM
                    );

                    if (min_distance > current_distance)
                    {
                        min_distance = current_distance;
                        min_cluster = clusters[j];
                    }
                }

                if (closest_centers[i] != min_cluster)
                {
                    // We changed a cluster-point pair
                    differ = 1;
                    closest_centers[i] = min_cluster;
                }

                // Update buffer that will be used to recalculate center
                //omp_set_lock(&lock_cluster_buffer);
                closest_centers[i]->add_to_buffer(points[i]->get_coordinates(), DIM);
                //omp_unset_lock(&lock_cluster_buffer);
            }
        }

        // If we found a difference from previous calculation
        if (!differ)
            break;

        // Get center for each cluster
        for (int i = 0; i < CLUSTERS; i++)
        {
            if (clusters[i]->get_points() == 0)
            {
                continue;
            }
            for (auto j = 0; j < DIM; j++)
            {
                clusters[i]->get_coordinates()[j] = clusters[i]->get_center_buffer()[j]
                                                    / clusters[i]->get_points();
            }
        }
    } while (differ && iterations < MAXITER);

    cout << endl << "Result: " << endl;
    cout << clusters[0]->get_coordinates()[0] << "  ";
    cout << clusters[0]->get_coordinates()[1];
}

int main()
{
	// Initialize seed
	srand((unsigned)time(0));

	// Generate points
	Point* points[POINTS];

    // Define some locks
    omp_lock_t cout_lock;
    omp_init_lock(&cout_lock);

    #pragma omp parallel num_threads(10)
    {
        #pragma omp for
        for (int i = 0; i < POINTS; i++)
        {
            //omp_set_lock(&cout_lock);
            //cout << "This is point " << i << endl;
            //omp_unset_lock(&cout_lock);

            double *coordinates = new double[DIM];
            for (int j = 0; j < DIM; j++)
            {
                coordinates[j] = (rand()%100)+1;
            }
            points[i] = new Point(coordinates);
        }
    }

	// Create cluster centers
	Cluster* clusters[CLUSTERS];
    Cluster* clusters_cpy[CLUSTERS];

    #pragma omp parallel num_threads(5)
    {
        #pragma omp for
        for (auto i = 0; i < CLUSTERS; i++)
        {
            // Create coordinates array and also a copy for multithreaded use
            double *coordinates = new double[DIM];
            double *coordinates_cpy = new double[DIM];

            for (auto j = 0; j < DIM; j++)
            {
                coordinates[j] = (rand() % 100) + 1;
                coordinates_cpy[j] = coordinates[j];
            }
            clusters[i] = new Cluster(coordinates);
            clusters_cpy[i] = new Cluster(coordinates_cpy);
        }
    }

    // Save execution times for before and after here
    high_resolution_clock::time_point t1,t2;
    long long int duration;

    // Single threaded run
    t1 = high_resolution_clock::now();
    single_threaded(clusters, points);
    t2 = high_resolution_clock::now();
    duration = duration_cast<milliseconds>( t2 - t1 ).count();
    cout << endl << "Single Threaded Execution duration : " << duration << endl << endl;

    // Multi threaded run
    t1 = high_resolution_clock::now();
    multi_threaded(clusters_cpy, points);
    t2 = high_resolution_clock::now();
    duration = duration_cast<milliseconds>( t2 - t1 ).count();
    cout << endl << "Multi Threaded Execution duration : " << duration;

    return 0;
}