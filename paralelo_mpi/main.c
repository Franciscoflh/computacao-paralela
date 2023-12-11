#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <mpi.h>

#define NUM_ANTS 50
#define NUM_ITERATIONS 100
#define PHEROMONE_EVAPORATION 0.5
#define ALPHA 1.0
#define BETA 2.0
#define NUM_CITIES 100

double pheromone_matrix[NUM_CITIES][NUM_CITIES];
double distance_matrix[NUM_CITIES][NUM_CITIES];

void calculate_probabilities(int current_city, int allowed_cities[], double probabilities[], double pheromone_matrix[NUM_CITIES][NUM_CITIES], double distance_matrix[NUM_CITIES][NUM_CITIES]) {
  double total_prob = 0.0;
  for (int i = 0; i < NUM_CITIES; i++) {
    if (allowed_cities[i] != -1) {
      double pheromone = pheromone_matrix[current_city][allowed_cities[i]];
      double distance = distance_matrix[current_city][allowed_cities[i]];
      probabilities[i] = pow(pheromone, ALPHA) * pow(1.0 / distance, BETA);
      total_prob += probabilities[i];
    } else {
      probabilities[i] = 0.0;
    }
  }
  for (int i = 0; i < NUM_CITIES; i++) {
    if (allowed_cities[i] != -1) {
      probabilities[i] /= total_prob;
    }
  }
}

void ant_tour(int ant_path[], double pheromone_matrix[NUM_CITIES][NUM_CITIES], double distance_matrix[NUM_CITIES][NUM_CITIES]) {
  int current_city = rand() % NUM_CITIES;
  ant_path[0] = current_city;
  int allowed_cities[NUM_CITIES];
  for (int i = 0; i < NUM_CITIES; i++) {
    allowed_cities[i] = i;
  }
  allowed_cities[current_city] = -1;

  for (int i = 1; i < NUM_CITIES; i++) {
    double probabilities[NUM_CITIES];
    calculate_probabilities(current_city, allowed_cities, probabilities, pheromone_matrix, distance_matrix);

    double rand_value = (double)rand() / RAND_MAX;
    double prob_sum = 0.0;
    int next_city = -1;
    for (int j = 0; j < NUM_CITIES; j++) {
      if (allowed_cities[j] != -1) {
        prob_sum += probabilities[j];
        if (rand_value <= prob_sum) {
          next_city = allowed_cities[j];
          break;
        }
      }
    }

    ant_path[i] = next_city;
    allowed_cities[next_city] = -1;
    current_city = next_city;
  }
}

void update_pheromones(int ant_paths[NUM_ANTS][NUM_CITIES], double pheromone_matrix[NUM_CITIES][NUM_CITIES], double distance_matrix[NUM_CITIES][NUM_CITIES]) {
  for (int i = 0; i < NUM_CITIES; i++) {
    for (int j = 0; j < NUM_CITIES; j++) {
      pheromone_matrix[i][j] *= (1.0 - PHEROMONE_EVAPORATION);
    }
  }

  for (int k = 0; k < NUM_ANTS; k++) {
    for (int i = 0; i < NUM_CITIES - 1; i++) {
      int city_i = ant_paths[k][i];
      int city_j = ant_paths[k][i + 1];
      pheromone_matrix[city_i][city_j] += 1.0 / distance_matrix[city_i][city_j];
    }
    int last_city = ant_paths[k][NUM_CITIES - 1];
    int first_city = ant_paths[k][0];
    pheromone_matrix[last_city][first_city] += 1.0 / distance_matrix[last_city][first_city];
  }
}

int main(int argc, char **argv) {
  int rank, size;
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  printf("Número total de processos: %d\n", size);

  srand(time(NULL) + rank); // Seed único para cada processo

  if (rank == 0) {
    for (int i = 0; i < NUM_CITIES; i++) {
      for (int j = 0; j < NUM_CITIES; j++) {
        pheromone_matrix[i][j] = 1.0;
        distance_matrix[i][j] = hypot(i - j, i - j);
      }
    }
  }

  MPI_Bcast(&pheromone_matrix, NUM_CITIES * NUM_CITIES, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&distance_matrix, NUM_CITIES * NUM_CITIES, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  double start_time = MPI_Wtime();
  int local_num_ants = NUM_ANTS / size + (rank < NUM_ANTS % size);
  int ant_paths[NUM_ANTS][NUM_CITIES];

  for (int iteration = 0; iteration < NUM_ITERATIONS; iteration++) {
    for (int i = 0; i < local_num_ants; i++) {
      ant_tour(ant_paths[i + rank * local_num_ants], pheromone_matrix, distance_matrix);
    }

    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Allgather(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, ant_paths, NUM_CITIES * local_num_ants, MPI_INT, MPI_COMM_WORLD);

    if (rank == 0) {
      update_pheromones(ant_paths, pheromone_matrix, distance_matrix);
    }

    MPI_Bcast(&pheromone_matrix, NUM_CITIES * NUM_CITIES, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  }

  double end_time = MPI_Wtime();

  if (rank == 0) {
    double best_distance = __DBL_MAX__;
    int best_path[NUM_CITIES];

    // Encontrar o melhor caminho
    for (int i = 0; i < NUM_ANTS; i++) {
      double distance = 0.0;
      for (int j = 0; j < NUM_CITIES - 1; j++) {
        distance += distance_matrix[ant_paths[i][j]][ant_paths[i][j + 1]];
      }
      distance += distance_matrix[ant_paths[i][NUM_CITIES - 1]][ant_paths[i][0]]; // Volta à cidade inicial

      if (distance < best_distance) {
        best_distance = distance;
        for (int k = 0; k < NUM_CITIES; k++) {
          best_path[k] = ant_paths[i][k];
        }
      }
    }

    printf("Melhor rota encontrada:\n");
    for (int i = 0; i < NUM_CITIES; i++) {
      printf("%d ", best_path[i]);
    }
    printf("\nDistância total: %f\n", best_distance);
    printf("Tempo de execução: %f segundos\n", end_time - start_time);
  }

  MPI_Finalize();
  return 0;
}
