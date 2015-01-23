import subprocess

# grid.tune("dof-dragons.dgauss.pbrt", [16], [128], [10], [0.0001, 0.01, 1.0, 100.0, 10000.0], [1.0])

def tune(input_file, samp_ress, reso_ress, search_regs, c1s, c2s):
    with open(input_file, "r") as f:
        inp = f.read()
    for samp_res in samp_ress:
        for reso_res in reso_ress:
            for search_reg in search_regs:
                for c1 in c1s:
                    for c2 in c2s:
                        out = inp
                        out = out.replace("$SAMPLE_RESO", str(samp_res))
                        out = out.replace("$RECONS_RESO", str(reso_res))
                        out = out.replace("$SEARCH_RANGE", str(search_reg))
                        out = out.replace("$C1", str(c1))
                        out = out.replace("$C2", str(c2))
                        output_file = "%s-%d-%d-%d-%s-%s.pbrt" % (input_file[:-5], samp_res, reso_res, search_reg, str(c1), str(c2))
                        with open(output_file, "w") as fout:
                            fout.write(out)
                        print "Running", output_file
                        subprocess.check_call(["./bin/pbrt.exe", output_file])
                        subprocess.check_call(["./exrtopng", output_file[:-5] + ".exr", output_file[:-5] + ".png"])
