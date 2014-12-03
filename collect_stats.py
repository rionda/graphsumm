import math
import os.path
import sys

def lastline(filename):
    if os.path.getsize(filename) >= 2:
        with open(filename, "rb") as f:
            f.seek(-2, 2)            # Jump to the second last byte.
            while f.read(1) != b'\n': # Until EOL is found...
                f.seek(-2, 1)        # ...jump back the read byte plus one more.
            last = f.readline()      # Read last line.
        return last.decode("ascii")
    else:
        sys.stderr.write("Empty file {}\n".format(filename))
        return ""

def main():
    float_stats_keywords = ["ind", "max_ind", "comp_size", "orig_size", "ratio",
            "L1_err", "L1_max_err", "L2_err", "L2_max_err", "cut-norm_err",
            "size_avg", "size_max", "size_min", "size_stddev", "dens_avg",
            "dens_max", "dens_min", "dens_stddev", "int_dens_avg",
            "int_dens_max", "int_dens_min", "int_dens_stddev", "ext_dens_avg",
            "ext_dens_max", "ext_dens_min", "ext_dens_stddev",
            "adj_err_avg", "adj_err_stdev", "adj_err_max", "adj_err_min",
            "deg_err_avg", "deg_err_stdev", "deg_err_max", "deg_err_min",
            "absdeg_err_avg", "absdeg_err_max", "absdeg_err_min",
            "reldeg_err_avg", "reldeg_err_stdev", "reldeg_err_max", "reldeg_err_min", 
            "clust_err" ]
    stats_keywords = ["alg", "in_file", "nodes", "edges", "avg_deg", "avg_dens",
            "triangles", "k", "approx_reconstr_err"] + float_stats_keywords
    summ_keywords = ["err_type", "dims", "dont_square", "random_summ",
    "mini_batch", "approx_alg", "total_time", "clustering_time"]
    grass_keywords = ["param_c"]

    # First argument is a list of comma separated values denoting which
    # statistics to print. If "-a", we print all of them.
    # Following arguments are input files.
    if len(sys.argv) < 3:
        sys.stderr.write("USAGE: {} [-a|list_of_csv_stats] file1 [file2 ...]\n".format(sys.argv[0]))
        sys.exit(1)
    if sys.argv[1] != "-a":
        stats_to_print = sys.argv[1].split(",")
    else:
        stats_to_print = []

    for stat in stats_to_print:
        if stat not in stats_keywords and stat not in summ_keywords and stat not in grass_keywords:
            sys.stderr.write("ERROR: keyword {} not recognized\n".format(stat))
            sys.exit(1)

    stats_sums = dict(zip(float_stats_keywords, [0] * len(float_stats_keywords)))
    squared_stats_sums = dict(zip(float_stats_keywords, [0] * len(float_stats_keywords)))
    stats_maxs = dict(zip(float_stats_keywords, [-10000000] * len(float_stats_keywords)))
    stats_mins = dict(zip(float_stats_keywords, [10000000] * len(float_stats_keywords)))

    files_num = 0
    stats = []
    for filename in sys.argv[2:]:
        stats_line = lastline(filename)
        if len(stats_line) > 1:
            files_num += 1
        else:
            sys.stderr.write("Empy line, skipping\n")
            continue
        curr_stats_list = stats_line.split(", ")
        # Fill key list with placeholders for additional stats
        curr_stats_keywords = stats_keywords.copy()
        for i in range(len(curr_stats_list) - len(stats_keywords)):
            curr_stats_keywords.append("addit_" + str(i))

        curr_stats = dict(zip(curr_stats_keywords, curr_stats_list))

        if curr_stats["alg"] == "summ":
            curr_stats["err_type"] = curr_stats.pop("addit_0")
            curr_stats["dims"] = curr_stats.pop("addit_1")
            curr_stats["dont_square"] = curr_stats.pop("addit_2")
            curr_stats["random_summ"] = curr_stats.pop("addit_3")
            curr_stats["mini_batch"] = curr_stats.pop("addit_4")
            curr_stats["approx_alg"] = curr_stats.pop("addit_5")
            curr_stats["total_time"] = float(curr_stats.pop("addit_6"))
            curr_stats["clustering_time"] = float(curr_stats.pop("addit_7"))
            curr_stats["squaring_time"] = float(curr_stats.pop("addit_8"))
            if "total_time" not in float_stats_keywords:
                for add_stat in ["total_time", "clustering_time", "squaring_time"]:
                    float_stats_keywords.append(add_stat)
                    stats_sums[add_stat] = 0
                    squared_stats_sums[add_stat] = 0
                    stats_maxs[add_stat] = -1000000
                    stats_mins[add_stat] = 1000000
        elif curr_stats["alg"] == "grass":
            curr_stats["err_type"] = curr_stats.pop("addit_0")
            curr_stats["param_c"] = curr_stats.pop("addit_1")
            curr_stats["total_time"] = float(curr_stats.pop("addit_2"))
            if "total_time" not in float_stats_keywords:
                float_stats_keywords += ["total_time"]
                stats_sums["total_time"] = 0
                squared_stats_sums["total_time"] = 0
                stats_maxs["total_time"] = -10000000 
                stats_mins["total_time"] = 10000000
        else:
            curr_stats["total_time"] = float(curr_stats.pop("addit_0"))
            if "total_time" not in float_stats_keywords:
                float_stats_keywords += ["total_time"]
                stats_sums["total_time"] = 0
                squared_stats_sums["total_time"] = 0
                stats_maxs["total_time"] = -10000000
                stats_mins["total_time"] = 10000000 

        curr_stats["L2_err"] = math.sqrt(float(curr_stats["L2_err"]))

        for stat in float_stats_keywords:
            curr_stats[stat] = float(curr_stats[stat])
            curr = curr_stats[stat]
            stats_sums[stat] += curr
            squared_stats_sums[stat] += math.pow(curr, 2.0)
            stats_maxs[stat] = max(stats_maxs[stat], curr)
            stats_mins[stat] = min(stats_mins[stat], curr)
        stats.append(curr_stats)

    if len(stats) == 0:
        sys.stderr.write("All files were empty. Printing fake lines and exiting.\n")
        print("")
        print("")
        sys.exit(0)

    if sys.argv[1] == "-a":
        stats_to_print = stats[0].keys()

    headings_to_print = []
    to_print = []
    for stat in stats_to_print:
        if "avg" in stat:
            if stat == "avg_dens" or stat == "avg_deg":
                to_print.append(stats[0][stat])
            else: 
                to_print.append(str(stats_sums[stat] / files_num))
            headings_to_print.append(stat)
        elif "stddev" in stat or "stdev" in stat: # THIS IS FUN
            if stat == "degree_err_stdev": 
                divider = int(stats[0]["nodes"])
            elif stat == "adj_err_stdev":
                divider = (int(stats[0]["nodes"]) * (int(stats[0]["nodes"]) -1)) / 2
            elif stat == "size_stddev" or stat == "int_dens_stddev":
                divider = int(stats[0]["k"])
            elif stat == "dens_stddev":
                divider = math.pow(int(stats[0]["k"]), 2) 
            elif stat == "ext_dens_stddev":
                divider = int(stats[0]["k"]) * (int(stats[0]["k"]) -1)
            squared_sum = 0;
            for curr_stats in stats:
                curr_stddev = curr_stats[stat]
                curr_var = math.pow(curr_stddev, 2)
                curr_squared_sum = (curr_var +
                        math.pow(curr_stats[stat.replace("stddev", "avg")], 2)) * divider
                squared_sum += curr_squared_sum
            variance = (squared_sum / (divider * files_num)) -  \
                math.pow(stats_sums[stat.replace("stddev", "avg")] / files_num, 2)
            if variance < 0 and variance > -1e-15: # assuming precision error
                variance = 0
            stddev = math.sqrt(variance)
            to_print.append(str(stddev))
            headings_to_print.append(stat)
        elif "max" in stat:
            to_print.append(str(stats_maxs[stat]))
            headings_to_print.append(stat)
        elif "min" in stat:
            if stat == "mini_batch":
                to_print.append(stats[0][stat])
            else:
                to_print.append(str(stats_mins[stat]))
            headings_to_print.append(stat)
        elif ("err" in stat and stat != "err_type" and stat != "approx_reconstr_err" and "degree" not in stat) or "time" in stat:
            to_print.append(str(stats_sums[stat] / files_num))
            headings_to_print.append(stat + "_avg")
            variance = (squared_stats_sums[stat] / files_num) - \
                math.pow(stats_sums[stat] / files_num, 2)
            if variance < 0 and variance > -1e-4: # assuming precision error
                variance = 0
            to_print.append(str(math.sqrt(variance)))
            headings_to_print.append(stat + "_stddev")
            to_print.append(str(stats_maxs[stat]))
            headings_to_print.append(stat + "_max")
            to_print.append(str(stats_mins[stat]))
            headings_to_print.append(stat + "_min")
        elif stat == "in_file":
            to_print.append(os.path.basename(stats[0][stat]).split(".")[0])
            headings_to_print.append(stat)
        elif stat == "approx_alg":
            if stats[0][stat] == "1":
                to_print.append("appr")
            else:
                to_print.append("gree")
            headings_to_print.append("algo")
        else:
            to_print.append(str(stats[0][stat]))
            headings_to_print.append(stat)

    print(",".join(headings_to_print))
    print(",".join(to_print))


if __name__ == "__main__":
    main()

