#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define NUM_ANTS 50
#define NUM_CITIES 500
#define NUM_ITERATIONS 100
#define PHEROMONE_EVAPORATION 0.5
#define ALPHA 1.0
#define BETA 2.0

int distance_matrix[NUM_CITIES][NUM_CITIES];
double pheromone_matrix[NUM_CITIES][NUM_CITIES];
int **ant_paths; 

int distance_matrix[NUM_CITIES][NUM_CITIES];
double pheromone_matrix[NUM_CITIES][NUM_CITIES];

void calculate_probabilities(int current_city, int allowed_cities[], double probabilities[]) {
    double total_prob = 0.0;
    for (int i = 0; i < NUM_CITIES; i++) {
        if (allowed_cities[i] != -1) {
            double pheromone = pheromone_matrix[current_city][allowed_cities[i]];
            double distance = distance_matrix[current_city][allowed_cities[i]];
            probabilities[i] = pow(pheromone, ALPHA) * pow(1.0 / distance, BETA);
            total_prob += probabilities[i];
        }
    }
    for (int i = 0; i < NUM_CITIES; i++) {
        if (allowed_cities[i] != -1) {
            probabilities[i] /= total_prob;
        }
    }
}

void ant_tour(int ant_path[]) {
    int current_city = rand() % NUM_CITIES;
    ant_path[0] = current_city;
    int allowed_cities[NUM_CITIES];
    for (int i = 0; i < NUM_CITIES; i++) {
        allowed_cities[i] = i;
    }
    allowed_cities[current_city] = -1;

    for (int i = 1; i < NUM_CITIES; i++) {
        double probabilities[NUM_CITIES] = {0.0};
        calculate_probabilities(current_city, allowed_cities, probabilities);

        double rand_value = (double)rand() / RAND_MAX;
        double prob_sum = 0.0;
        int next_city = -1;
        for (int j = 0; j < NUM_CITIES && next_city == -1; j++) {
            if (allowed_cities[j] != -1) {
                prob_sum += probabilities[j];
                if (rand_value <= prob_sum) {
                    next_city = allowed_cities[j];
                }
            }
        }

        if (next_city == -1) {
            next_city = allowed_cities[0];
        }

        ant_path[i] = next_city;
        allowed_cities[next_city] = -1;
        current_city = next_city;
    }
}

void update_pheromones(int **ant_paths) {
    double delta_pheromones[NUM_CITIES][NUM_CITIES] = {0};

    for (int i = 0; i < NUM_ANTS; i++) {
        for (int j = 0; j < NUM_CITIES - 1; j++) {
            delta_pheromones[ant_paths[i][j]][ant_paths[i][j + 1]] += 1.0 / distance_matrix[ant_paths[i][j]][ant_paths[i][j + 1]];
        }
        delta_pheromones[ant_paths[i][NUM_CITIES - 1]][ant_paths[i][0]] += 1.0 / distance_matrix[ant_paths[i][NUM_CITIES - 1]][ant_paths[i][0]];
    }

    for (int i = 0; i < NUM_CITIES; i++) {
        for (int j = 0; j < NUM_CITIES; j++) {
            pheromone_matrix[i][j] *= (1.0 - PHEROMONE_EVAPORATION);
            pheromone_matrix[i][j] += delta_pheromones[i][j];
        }
    }
}

int main() {
    srand(time(NULL));

    // Alocando memória dinamicamente para ant_paths
    ant_paths = (int **)malloc(NUM_ANTS * sizeof(int *));
    for (int i = 0; i < NUM_ANTS; i++) {
        ant_paths[i] = (int *)malloc(NUM_CITIES * sizeof(int));
    }

    // Inicialização das matrizes de distância e feromônio
    for (int i = 0; i < NUM_CITIES; i++) {
        for (int j = 0; j < NUM_CITIES; j++) {
            if (i == j) {
                distance_matrix[i][j] = 0;
            } else {
                distance_matrix[i][j] = rand() % 100 + 1;  // Distâncias aleatórias entre 1 e 100
            }
            pheromone_matrix[i][j] = 1.0;
        }
    }

    clock_t start_time = clock();

    for (int iteration = 0; iteration < NUM_ITERATIONS; iteration++) {
        for (int i = 0; i < NUM_ANTS; i++) {
            ant_tour(ant_paths[i]);
        }

        update_pheromones(ant_paths);
    }

    int best_path[NUM_CITIES];
    double best_distance = __DBL_MAX__;

    for (int i = 0; i < NUM_ANTS; i++) {
        double distance = 0.0;
        for (int j = 0; j < NUM_CITIES - 1; j++) {
            distance += distance_matrix[ant_paths[i][j]][ant_paths[i][j + 1]];
        }
        distance += distance_matrix[ant_paths[i][NUM_CITIES - 1]][ant_paths[i][0]];

        if (distance < best_distance) {
            best_distance = distance;
            for (int j = 0; j < NUM_CITIES; j++) {
                best_path[j] = ant_paths[i][j];
            }
        }
    }

    clock_t end_time = clock();
    double time_taken = (double)(end_time - start_time) / CLOCKS_PER_SEC;

    printf("Melhor rota encontrada: ");
    for (int i = 0; i < 10 && i < NUM_CITIES; i++) {  // Imprimindo apenas as primeiras 10 cidades para brevidade
        printf("%d ", best_path[i]);
    }
    printf("...");
    printf("\nDistância total: %.2f\n", best_distance);
    printf("Tempo de execução: %.2f segundos\n", time_taken);

    // Liberando a memória alocada dinamicamente
    for (int i = 0; i < NUM_ANTS; i++) {
        free(ant_paths[i]);
    }
    free(ant_paths);

    return 0;
}