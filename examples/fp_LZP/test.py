def read_outfp(fname='out.fp'):
    with open(fname,'r') as f:
        lines = f.readlines()
    return lines

out = read_outfp('out.fp')
outref = read_outfp('out.fp.REF')

#...Num of lines changed since 2021-06-30
# assert len(out) == len(outref)

Ls = []
for l in out:
    assert 'nan' not in l or 'NaN' not in l

    dat = l.split()
    if 'iid' in l:
        Ls.append(float(dat[6]))
    if 'step,time,best,vars' in l:
        dat = l.split()
        best = float(dat[3])
        assert abs(best -min(Ls)) < 0.001

print('pass')
