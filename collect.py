import sys
import math

def main():
    if len(sys.argv) != 2:
        sys.stderr.write("USAGE: {} file\n".format(sys.argv[0]))
        sys.exit(1)

    with open(sys.argv[1], 'rt') as resfile:
        for line in resfile:
            lines = [line.split(",")]
            for i in range(4):
                lines.append(resfile.readline().split(","))
            adj_sum = 0
            adj_stdev_sum = 0
            adj_max = -1
            adj_min = 2
            deg_sum = 0
            deg_stdev_sum = 0
            deg_max = -1
            deg_min = 2
            clust_sum = 0
            clust_stdev_sum = 0
            clust_max = -1
            clust_min = 2
            for i in range(5):
                adj_sum += float(lines[i][3])
                adj_stdev_sum += float(lines[i][4])
                adj_max = max(adj_max, float(lines[i][5]))
                adj_min = min(adj_min, float(lines[i][6]))
                deg_sum += float(lines[i][7])
                deg_stdev_sum += float(lines[i][8])
                deg_max = max(deg_max, float(lines[i][9]))
                deg_min = min(deg_min, float(lines[i][10]))
                clust_sum += float(lines[i][11])
                clust_stdev_sum += pow(float(lines[i][11]), 2);
                clust_max = max(clust_max, float(lines[i][11]))
                clust_min = min(clust_min, float(lines[i][11]))
            adj_avg = adj_sum / 5
            adj_stdev_avg = adj_stdev_sum / 5
            deg_avg = deg_sum / 5
            deg_stdev_avg = deg_stdev_sum / 5
            clust_avg = clust_sum / 5
            clust_stdev = math.sqrt(clust_stdev_sum / 5 - pow(clust_avg, 2))
            print("{}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}".format(
                lines[0][0], lines[0][1], lines[0][2], adj_avg, adj_stdev_avg,
                adj_max, adj_min, deg_avg, deg_stdev_avg, deg_max, deg_min,
                clust_avg, clust_stdev, clust_max, clust_min))
    return 0

if __name__ == "__main__":
    main()

