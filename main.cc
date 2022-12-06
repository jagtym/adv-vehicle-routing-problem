#include <algorithm>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <random>
#include <sstream>
#include <vector>

using namespace std;

#include "customer.h"
#include "truck.h"

#define DEPOT 0
#define MIN_TIME 99999999

random_device dev;
mt19937 rng(dev());

int random(int lower_bound, int higher_bound) {
    uniform_int_distribution<mt19937::result_type> dist(lower_bound, higher_bound);
    return dist(rng);
}

double distance_sum(vector<Truck> solution) {
    double distance = 0;
    for (int i = 0; i < solution.size(); i++) {
        distance += solution[i].act_time;
    }
    return distance;
}

void load_file_contents(ifstream &file, vector<string> &file_content) {
    string line;
    while (file >> ws && getline(file, line)) {
        file_content.push_back(line);
    }
}

void write_to_file(ofstream &output_file, vector<Truck> &trucks) {
    if (trucks.size() == 0) {
        output_file << -1 << " " << 0 << endl;
    } else {
        double distance = 0;
        for (int i = 0; i < trucks.size(); i++) {
            distance += trucks[i].act_time;
        }
        cout << "WRITE TO FILE: " << trucks.size() << " " << fixed << setprecision(5) << distance << endl;
        // cout << "Distance sum: " << distance_sum(trucks) << endl;
        output_file << trucks.size() << " " << fixed << setprecision(5) << distance << endl;
        for (int i = 0; i < trucks.size(); i++) {
            for (int j = 0; j < trucks[i].customers.size(); j++) {
                output_file << trucks[i].customers[j] << " ";
            }
            output_file << endl;
        }
    }
}

float get_distance(int x1, int y1, int x2, int y2) {
    return sqrt((float)pow(x1 - x2, 2) + (float)pow(y1 - y2, 2));
}

void populate_customer_array(vector<string> &file_content, vector<Customer> &customers) {
    const int first_customer_line = 6;

    for (int i = first_customer_line; i < file_content.size(); i++) {
        int index, x_cord, y_cord, demand, ready_time, due_date, service_time;

        istringstream(file_content[i]) >> index >> x_cord >> y_cord >> demand >> ready_time >> due_date >> service_time;
        Customer customer = Customer(index, x_cord, y_cord, demand, ready_time, due_date, service_time);
        customers.push_back(customer);
    }
}

void populate_matrix(vector<Customer> &customers, vector<vector<float>> &matrix) {
    for (int y = 0; y < customers.size(); y++) {
        for (int x = 0; x < customers.size(); x++) {
            int first_customer_x = customers[x].x_cord;
            int first_customer_y = customers[x].y_cord;
            int second_customer_x = customers[y].x_cord;
            int second_customer_y = customers[y].y_cord;
            matrix[y][x] = get_distance(first_customer_x, first_customer_y, second_customer_x, second_customer_y);
        }
    }
}

void find_answer(vector<Customer> &customers, vector<vector<float>> &matrix, vector<Truck> &trucks, const int capacity) {
    int minimum_time = MIN_TIME;
    int visited_customers = 0;
    int next_location = DEPOT;
    customers[0].is_served = true;

    Truck first_truck = Truck(capacity);
    trucks.push_back(first_truck);

    while (visited_customers < customers.size() - 1) {
        int current_truck = trucks.size() - 1;
        for (int i = 0; i < customers.size(); i++) {
            float customer_distance = matrix[trucks[current_truck].act_location][i];
            float depot_distance = matrix[i][0];
            Customer depot = customers[0];
            Customer current_customer = customers[i];

            if (i == trucks[current_truck].act_location) {
                continue;
            }
            if (current_customer.is_served == true) {
                continue;
            }
            if (trucks[current_truck].act_capacity < current_customer.demand) {
                continue;
            }
            if (trucks[current_truck].act_time + customer_distance > current_customer.due_date) {
                continue;
            }
            if (trucks[current_truck].act_time + customer_distance + current_customer.service_time + depot_distance > depot.due_date) {
                continue;
            }
            if (current_customer.ready_time + current_customer.service_time + depot_distance > depot.due_date) {
                continue;
            }

            if (customers[i].ready_time + customers[i].service_time < minimum_time) {
                minimum_time = customers[i].ready_time + customers[i].service_time;
                next_location = i;
            }
        }

        if (next_location != DEPOT) {
            trucks[current_truck].customers.push_back(next_location);
            Customer current_customer = customers[next_location];
            float next_distance = matrix[trucks[current_truck].act_location][next_location];

            if (next_distance + trucks[current_truck].act_time < current_customer.ready_time) {
                trucks[current_truck].act_time += current_customer.ready_time + current_customer.service_time - trucks[current_truck].act_time;
                trucks[current_truck].act_capacity -= current_customer.demand;
                trucks[current_truck].act_location = next_location;
            } else {
                trucks[current_truck].act_time += next_distance + current_customer.service_time;
                trucks[current_truck].act_capacity -= current_customer.demand;
                trucks[current_truck].act_location = next_location;
            }

            visited_customers++;
            customers[next_location].is_served = true;
            next_location = DEPOT;
            minimum_time = MIN_TIME;
        } else if (trucks[current_truck].customers.size() == 0) {
            trucks.clear();
            break;
        } else {
            float depot_distance = matrix[trucks[current_truck].act_location][0];
            trucks[current_truck].act_time += depot_distance;
            minimum_time = MIN_TIME;
            next_location = DEPOT;
            Truck new_truck = Truck(capacity);
            trucks.push_back(new_truck);
        }
    }
    trucks[trucks.size() - 1].act_time += matrix[trucks[trucks.size() - 1].act_location][0];
}

void sort_neighbourhood(vector<vector<Truck>> &neighbourhood) {
    vector<vector<Truck>> prev_neighbourhood;
    for (int i = 0; i < neighbourhood.size(); i++) {
        prev_neighbourhood.push_back(neighbourhood[i]);
    }
    neighbourhood.clear();
    for (int i = 0; i < prev_neighbourhood.size(); i++) {
        int min_trucks = 9999999;
        float min_distance = 9999999;
        int index = 0;
        for (int j = 0; j < prev_neighbourhood[i].size(); j++) {
            if (prev_neighbourhood[i].size() < min_trucks) {
                min_trucks = prev_neighbourhood[i].size();
                min_distance = distance_sum(prev_neighbourhood[i]);
                index = i;
            } else if (prev_neighbourhood.size() == min_trucks && distance_sum(prev_neighbourhood[i]) < min_distance) {
                min_distance = distance_sum(prev_neighbourhood[i]);
                index = i;
            }
        }
        neighbourhood.push_back(prev_neighbourhood[index]);
        prev_neighbourhood.erase(prev_neighbourhood.begin() + index);
    }
    return;
}

bool compare_solutions(vector<Truck> t1, vector<Truck> t2) {
    if (t1.size() > t2.size()) {
        return true;
    }
    if (t1.size() == t2.size() && distance_sum(t1) > distance_sum(t2)) {
        return true;
    }
    return false;
}

bool compare_solutions2(vector<Truck> &a, vector<Truck> &b) {
    if (a.size() != b.size())
        return a.size() < b.size();
    return distance_sum(a) < distance_sum(b);
}

void get_neighbourhood(vector<Customer> &customers, vector<vector<float>> &matrix, vector<Truck> &current_solution, vector<vector<Truck>> &neighbourhood, int capacity) {
    int vehicle_nr = random(0, current_solution.size() - 1);
    int customer_nr_on_route = random(0, current_solution[vehicle_nr].customers.size() - 1);
    int customer_nr = current_solution[vehicle_nr].customers[customer_nr_on_route];
    for (int i = 0; i < current_solution.size(); i++) {
        if (i == vehicle_nr) {
            continue;
        }
        if (customers[customer_nr].demand > current_solution[i].act_capacity) {
            continue;
        }

        int real_time, possible_time = 0;
        int current_position = 0;
        vector<Truck> possible_solution;
        Truck possible_new_truck(capacity);
        bool pushed = false;
        for (int j = 0; j < current_solution[i].customers.size(); j++) {
            if (!pushed) {
                real_time += max(matrix[current_position][current_solution[i].customers[j]], (float)customers[current_solution[i].customers[j]].ready_time - real_time) + customers[current_solution[i].customers[j]].service_time;
                possible_time += max(matrix[current_position][customer_nr], (float)customers[customer_nr].ready_time - possible_time) + customers[customer_nr].service_time;
                possible_time += max(matrix[customer_nr][j], (float)customers[current_solution[i].customers[j]].ready_time - possible_time) + customers[current_solution[i].customers[j]].service_time;
                if (possible_time <= real_time) {
                    possible_new_truck.customers.push_back(customer_nr);
                    pushed = true;
                } else {
                    possible_time = real_time;
                }
                current_position = current_solution[i].customers[j];
            }
            possible_new_truck.customers.push_back(current_solution[i].customers[j]);
        }

        if (pushed) {
            for (int j = 0; j < current_solution.size(); j++) {
                if (i != j) {
                    possible_solution.push_back(current_solution[j]);
                } else {
                    possible_solution.push_back(possible_new_truck);
                }
            }
            possible_solution[vehicle_nr].customers.erase(possible_solution[vehicle_nr].customers.begin() + customer_nr_on_route);

            // changing act_time in line where customer was deleted and where customer was put and act_capacity
            bool validate_route = true;

            possible_solution[vehicle_nr].act_time = 0;
            possible_solution[vehicle_nr].act_capacity = capacity;
            current_position = 0;
            for (int j = 0; j < possible_solution[vehicle_nr].customers.size(); j++) {
                possible_solution[vehicle_nr].act_time += max(matrix[current_position][possible_solution[vehicle_nr].customers[j]], customers[possible_solution[vehicle_nr].customers[j]].ready_time - possible_solution[vehicle_nr].act_time) + customers[possible_solution[vehicle_nr].customers[j]].service_time;
                possible_solution[vehicle_nr].act_capacity -= customers[possible_solution[vehicle_nr].customers[j]].demand;
                current_position = possible_solution[vehicle_nr].customers[j];
            }
            possible_solution[vehicle_nr].act_time += matrix[current_position][0];

            possible_solution[i].act_time = 0;
            possible_solution[i].act_capacity = capacity;
            current_position = 0;
            for (int j = 0; j < possible_solution[i].customers.size(); j++) {
                possible_solution[i].act_time += max(matrix[current_position][possible_solution[i].customers[j]], customers[possible_solution[i].customers[j]].ready_time - possible_solution[i].act_time) + customers[possible_solution[i].customers[j]].service_time;
                possible_solution[i].act_capacity -= customers[possible_solution[i].customers[j]].demand;
                if (possible_solution[i].act_time - customers[possible_solution[i].customers[j]].service_time > customers[possible_solution[i].customers[j]].due_date) {
                    validate_route = false;
                }
                current_position = possible_solution[i].customers[j];
            }
            possible_solution[i].act_time += matrix[current_position][0];
            if (possible_solution[i].act_time > customers[0].due_date) {
                validate_route = false;
            }

            if (possible_solution[vehicle_nr].customers.size() == 0) {
                // cout << "Empty truck" << endl;
                possible_solution.erase(possible_solution.begin() + vehicle_nr);
            }

            if (validate_route) {
                neighbourhood.push_back(possible_solution);
            }
        }
    }
    return;
}

bool in_tabu2(vector<vector<Truck>> &tabu_list, vector<Truck> &solution) {
    auto compare_trucks = [](Truck t1, Truck t2) {
        return equal(t1.customers.begin(), t1.customers.end(), t2.customers.begin());
    };

    for (vector<Truck> tabu : tabu_list) {
        if (tabu.size() != solution.size()) {
            continue;
        }

        if (equal(tabu.begin(), tabu.end(), solution.begin(), compare_trucks)) {
            return true;
        }
    }
    return false;
}

bool in_tabu(vector<vector<Truck>> &tabu_list, vector<Truck> &solution) {
    for (int i = 0; i < tabu_list.size(); i++) {
        bool is_member = true;
        if (tabu_list[i].size() != solution.size()) {
            continue;
        }
        for (int x = 0; x < tabu_list[i].size(); x++) {
            if (!equal(tabu_list[i][x].customers.begin(), tabu_list[i][x].customers.end(), solution[x].customers.begin())) {
                is_member = false;
                break;
            }
        }
        if (is_member) {
            return true;
        }
    }
    return false;
}

bool new_solution_better(vector<Truck> &prev_solution, vector<Truck> &new_solution) {
    if (new_solution.size() < prev_solution.size()) {
        return true;
    }
    if (new_solution.size() == prev_solution.size() && distance_sum(new_solution) < distance_sum(prev_solution)) {
        return true;
    }
    return false;
}

void tabu_search(vector<Customer> &customers, vector<vector<float>> &matrix, vector<Truck> &trucks, int capacity, ofstream &output_file, int tabu_size) {
    vector<vector<Truck>> neigbourhood;
    vector<Truck> best_solution;
    vector<Truck> current_solution;
    for (int i = 0; i < trucks.size(); i++) {
        best_solution.push_back(trucks[i]);
        current_solution.push_back(trucks[i]);
    }

    vector<vector<Truck>> tabu_list;

    int how_many = 1;
    write_to_file(output_file, best_solution);

    while (true) {
        get_neighbourhood(customers, matrix, current_solution, neigbourhood, capacity);
        sort(neigbourhood.begin(), neigbourhood.end(), compare_solutions2);
        // cout << "Sort" << endl;
        // sort_neighbourhood(neigbourhood);
        // cout << "Sorted" << endl;
        vector<Truck> best_candidate;
        for (int i = 0; i < neigbourhood.size(); i++) {
            best_candidate.clear();
            if (!in_tabu2(tabu_list, neigbourhood[i])) {
                for (int j = 0; j < neigbourhood[i].size(); j++) {
                    best_candidate.push_back(neigbourhood[i][j]);
                }
                break;
            }
        }
        if (!best_candidate.empty()) {
            // cout << "Found candidate" << endl;
            tabu_list.push_back(best_candidate);
            if (tabu_list.size() > tabu_size) {
                tabu_list.erase(tabu_list.begin());
            }
            if (new_solution_better(best_solution, best_candidate)) {
                // cout << "New solution is better" << endl;
                best_solution.clear();
                current_solution.clear();
                for (int i = 0; i < best_candidate.size(); i++) {
                    best_solution.push_back(best_candidate[i]);
                    current_solution.push_back(best_candidate[i]);
                }
                output_file.close();
                output_file.open("tabu.txt");
                write_to_file(output_file, best_solution);
            } else {
                // Po zaczeciu brania gorszych ma czasem problem znalezc lepsze
                // cout << "New solution is worse" << endl;
                current_solution.clear();
                for (int i = 0; i < best_candidate.size(); i++) {
                    current_solution.push_back(best_candidate[i]);
                }
            }
            how_many--;
        }
        neigbourhood.clear();
    }
}

int main(int argc, char *argv[]) {
    auto start = std::chrono::high_resolution_clock::now();
    if (argc != 3) {
        cout << "usage: " << argv[0] << " dataset_file output_file" << endl;
        return 1;
    }
    const string input_filename = argv[1];
    const string output_filename = argv[2];

    ifstream input_file;
    ofstream output_file;
    input_file.open(input_filename);
    output_file.open(output_filename);

    if (!input_file.is_open() || !output_file.is_open()) {
        cerr << "File not found!" << endl;
        return 1;
    }

    string problem_name;
    int vehicle_number;
    int capacity;

    vector<string> file_content;
    load_file_contents(input_file, file_content);

    problem_name = file_content[0];
    istringstream(file_content[3]) >> vehicle_number >> capacity;

    vector<Customer> customers;
    populate_customer_array(file_content, customers);

    vector<vector<float>> matrix(customers.size(), vector<float>(customers.size(), 0));
    populate_matrix(customers, matrix);

    vector<Truck> trucks;

    find_answer(customers, matrix, trucks, capacity);
    write_to_file(output_file, trucks);

    int tabu_size = 10;

    ofstream output_tabu;
    output_tabu.open("tabu.txt");

    tabu_search(customers, matrix, trucks, capacity, output_tabu, tabu_size);
    if (output_tabu.is_open()) {
        output_tabu.close();
    }
    auto elapsed = std::chrono::high_resolution_clock::now() - start;
    long long microseconds = std::chrono::duration_cast<std::chrono::seconds>(elapsed).count();
    cout << "Czas: " << microseconds;
    input_file.close();
    output_file.close();
}