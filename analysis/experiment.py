import argparse
import csv
import matplotlib.pyplot as plt

def build_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument("experiment_output", type=argparse.FileType('r'))
    return parser;

def build_data(csv_reader):
    data = {} 
    print(csv_reader.fieldnames)
    for field in csv_reader.fieldnames:
        data[field] = []

    for row in csv_reader:
        for field in csv_reader.fieldnames:
            data[field].append(float(row[field]))
    return data

def plot_output(csv_file):
    data = build_data(csv.DictReader(csv_file))
    fig, axs = plt.subplots(nrows=2, ncols=4)
    i = 0;
    j = 0;
    for s in ["H", "P", "C"]:
        axs[i][j].plot(data["T"], data[s], label=s + " (mol/l)")
        axs[i][j].legend()
        j = j+1
    i = 1
    j = 0
    for s in ["S", "SH", "SP", "SC"]:
        axs[i][j].plot(data["T"], data[s], label=s + "/N")
        axs[i][j].legend()
        j = j+1
    plt.show()

if __name__ == "__main__":
    parser = build_parser();
    args = parser.parse_args()

    plot_output(args.experiment_output)

