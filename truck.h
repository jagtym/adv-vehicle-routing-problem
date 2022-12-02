struct Truck {
    int capacity;
    int act_capacity;
    int act_location;
    int act_time;
    vector<int> customers;

    Truck(int capacity) : capacity(capacity),
                          act_capacity(capacity),
                          act_location(0),
                          act_time(0) {}

    void visit(int act_time, int act_capacity, int act_location) {
        this->act_capacity -= act_capacity;
        this->act_time += act_time;
        this->act_location = act_location;
    }
};