#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>

using namespace std;

#include "customer.h"
#include "truck.h"

void load_file_contents(ifstream &file, vector<string> &file_content) {
    string line;
    while (file >> ws && getline(file, line)) {
        file_content.push_back(line);
    }
}

float get_distance(int x1, int y1, int x2, int y2) {
    return sqrt((float)pow(x1 - x2, 2) + (float)pow(y1 - y2, 2));
}

void populate_customer_array(vector<string> &file_content, vector<Customer> &customers) {
    const int first_customer_line = 6;

    for (int i = first_customer_line; i < file_content.size(); i++) {
        int index;
        int x_cord;
        int y_cord;
        int demand;
        int ready_time;
        int due_date;
        int service_time;

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
    const int MIN_TIME = 99999999;

    int minimum_time = MIN_TIME;
    int visited_customers = 0;
    int next_location = 0;
    customers[0].is_served = true;

    Truck first_truck = Truck(capacity);
    trucks.push_back(first_truck);

    while (visited_customers < customers.size() - 1) {
        int current_truck = trucks.size() - 1;
        for (int i = 0; i < customers.size(); i++) {
            cout << i << endl;
            if (i == trucks[current_truck].act_location) {
                cout << "1 IF" << endl;
                continue;
            }
            if (customers[i].is_served == true) {
                cout << "2 IF" << endl;
                continue;
            }
            if (trucks[current_truck].act_capacity < customers[i].demand) {
                cout << trucks[current_truck].act_capacity << " " << customers[i].demand << endl;
                cout << "3 IF" << endl;
                continue;
            }
            if (trucks[current_truck].act_time + matrix[trucks[current_truck].act_location][i] > customers[i].due_date) {
                cout << "4 IF" << endl;
                continue;
            }
            if (trucks[current_truck].act_time + matrix[trucks[current_truck].act_location][i] + customers[i].service_time + matrix[i][0] > customers[0].due_date) {
                cout << "5 IF" << endl;
                continue;
            }
            if (customers[i].ready_time + customers[i].service_time + matrix[i][0] > customers[0].due_date) {
                cout << "6 IF" << endl;
                continue;
            }

            if (customers[i].ready_time + customers[i].service_time < minimum_time) {
                minimum_time = customers[i].ready_time + customers[i].service_time;
                next_location = i;
            }
        }
        if (next_location != 0) {
            trucks[current_truck].customers.push_back(next_location);
            // cout << "Porownanie: " << matrix[trucks[current_truck].act_location][next_location] << " < " << customers[next_location].ready_time << endl;
            if (matrix[trucks[current_truck].act_location][next_location] < customers[next_location].ready_time) {
                // cout << "Liczy: " << trucks[current_truck].act_time << " + " << customers[next_location].ready_time << "+" << customers[next_location].service_time << "-" << trucks[current_truck].act_time << endl;
                trucks[current_truck].visit(customers[next_location].ready_time + customers[next_location].service_time - trucks[current_truck].act_time, customers[next_location].demand, next_location);
            } else {
                trucks[current_truck].visit(matrix[trucks[current_truck].act_location][next_location] + customers[next_location].service_time, customers[next_location].demand, next_location);
            }
            visited_customers++;
            customers[next_location].is_served = true;
            next_location = 0;
        } else if (trucks[current_truck].customers.size() == 0) {
            trucks.clear();
            break;
        } else {
            // cout << trucks[current_truck].act_location << " " << trucks[current_truck].act_time << " " << matrix[trucks[current_truck].act_location][0] << endl;
            trucks[current_truck].act_time += matrix[trucks[current_truck].act_location][0];
            // cout << "Nr trucka: " << current_truck << " Act_time: " << trucks[current_truck].act_time << endl;
            minimum_time = MIN_TIME;
            next_location = 0;
            Truck new_truck = Truck(capacity);
            trucks.push_back(new_truck);
        }
    }
    trucks[trucks.size() - 1].act_time += matrix[trucks[trucks.size() - 1].act_location][0];
}

void write_to_file(ofstream &output_file, vector<Truck> &trucks) {
    if (trucks.size() == 0) {
        output_file << -1 << " " << 0 << endl;
    } else {
        double distance = 0;
        for (int i = 0; i < trucks.size(); i++) {
            distance += trucks[i].act_time;
        }
        output_file << trucks.size() << " " << fixed << setprecision(5) << distance << endl;
        for (int i = 0; i < trucks.size(); i++) {
            for (int j = 0; j < trucks[i].customers.size(); j++) {
                output_file << trucks[i].customers[j] << " ";
            }
            output_file << endl;
        }
    }
}

int main(int argc, char *argv[]) {
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

    input_file.close();
    output_file.close();
}