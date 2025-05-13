import sys
import matplotlib.pyplot as plt

trans = {
    "CH" : ["Eilice", "blue", '-'],
    "oCH" : ["Eilice (inexact)", "blue", ':'],
    "CQ" : ["Simplified OvL", "orange", '-'],
    "ED" : ["Eddy", "green", '-.'],
    "GA" : ["GA", "purple", '-.'],
    "JAVA": ["CHHN (inexact)", "red", ':'],
    "RK": ["Eilice (rank ordered)", "blue", '--']
}
def parse_file(file):
    res = {
        "name": file.name,
        "model": "",
        "Nums": [],
        "Inserting": [],
        "Removing": [],
    }
    acc = 0
    mode = ""
    for line in file:
        line = line.rstrip()
        if res["model"] == "":
            res["model"] = line
            print(line)
            continue

        if "Inserting" in line or "Removing" in line:
            mode = line
            acc = 0
            continue

        row = line.split(' ')
        if acc >= len(res["Nums"]):
            res["Nums"].append(int(row[0]))

        if acc >= len(res[mode]):
            if(res["model"] == "JAVA" or res["model"] == "GA"):
                res[mode].append([float(row[2])])
            else:
                res[mode].append([float(row[2])])
        else:
            if(res["model"] == "JAVA" or res["model"] == "GA"):
                res[mode][acc].append(float(row[2]))
            else:
                res[mode][acc].append(float(row[2]))
        acc += 1

    res["Inserting"] = [sorted(x)[len(x)//2] for x in res["Inserting"]]
    res["Removing"] = res["Removing"].sort()
    return res


def make_plot(data):
    plt.tight_layout()
    ls = ['-','--','-.',(0,(3,2,1,2,1,2)),':']
    cl = ["orange","blue","green"]
    fig, ax = plt.subplots()
    for d in data:
        ax.plot(d["Nums"],d["Inserting"],label = trans[d["model"]][0],c=trans[d["model"]][1],ls=trans[d["model"]][2])

    #ax.legend(loc="upper left")
    ax.set(xlabel="No. of points",ylabel="Time in seconds",title="Dynamic construction time on uniform data",xlim=[50000,1e6],ylim=[0,900])
    ax.legend()
    fig.savefig("construct_unif.png", dpi=300.)

def main():
    # Parse
    data = []
    for path in sys.stdin:
        file = open(path[:-1])
        data += [parse_file(file)]

    make_plot(data)


if __name__ == "__main__":
    main()
