struct Customer {
    int index;
    int x_cord;
    int y_cord;
    int demand;
    int ready_time;
    int due_date;
    int service_time;

    bool is_served;

    Customer(int i, int x, int y, int d, int rt, int dd, int st) : index(i),
                                                                   x_cord(x),
                                                                   y_cord(y),
                                                                   demand(d),
                                                                   ready_time(rt),
                                                                   due_date(dd),
                                                                   service_time(st),
                                                                   is_served(false) {
    }

    void print() {
        cout << index << " " << x_cord << " " << y_cord << " " << demand << " " << ready_time << " " << due_date << " " << service_time << endl;
    }
};