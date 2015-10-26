
f = open("run_template")
template = f.readlines()[1]
f.close()

for i in range(50):
    ti = template.replace("xxx", str(i+1))
    g = open("bash_scripts/run_" + str(i+1), "w")
    g.write(ti)
    g.close()

