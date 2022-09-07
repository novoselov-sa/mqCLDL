import sys

nrange = []
totalsqUgs = {}
totallgUgs = {}
totallgQgs = {}
totallgindex = {}
numsqUgs = {}
numlgUgs = {}
numlgQgs = {}
numlgindex = {}

while True:
  line = sys.stdin.readline()
  if not line: break
  x = line.strip().split()
  if x[0] == '%':
    if x[1] == 'Ugs':
      n = int(x[2])
      Ugs = RR(x[3])
      if not n in nrange: nrange += [n]
      if not n in totalsqUgs: totalsqUgs[n] = 0
      if not n in totallgUgs: totallgUgs[n] = 0
      if not n in numsqUgs: numsqUgs[n] = 0
      if not n in numlgUgs: numlgUgs[n] = 0
      lgUgs = log(Ugs)/log(2.0)
      totalsqUgs[n] += Ugs*Ugs; numsqUgs[n] += 1
      totallgUgs[n] += lgUgs; numlgUgs[n] += 1
    if x[1] == 'Qgs':
      n = int(x[2])
      if not n in nrange: nrange += [n]
      Qgs = RR(x[3])
      if not n in totallgQgs: totallgQgs[n] = 0
      if not n in numlgQgs: numlgQgs[n] = 0
      lgQgs = log(Qgs)/log(2.0)
      totallgQgs[n] += lgQgs; numlgQgs[n] += 1
    if x[1] == 'lgindex':
      n = int(x[2])
      if not n in nrange: nrange += [n]
      lgindex = RR(x[3])
      if not n in totallgindex: totallgindex[n] = 0
      if not n in numlgindex: numlgindex[n] = 0
      totallgindex[n] += lgindex; numlgindex[n] += 1

avgsqUgs = {}
avglgUgs = {}
avglgQgs = {}
avglgindex = {}
for n in nrange:
  avgsqUgs[n] = totalsqUgs[n] / numsqUgs[n]
  avglgUgs[n] = totallgUgs[n] / numlgUgs[n]
  avglgQgs[n] = totallgQgs[n] / numlgQgs[n]
  avglgindex[n] = totallgindex[n] / numlgindex[n]

nrange = sorted(nrange)
print('$n$&' + '&'.join('%d' % n for n in nrange) + '\\\\')
print('\\noalign{\\hrule}')
print('average $\\log_2||u^*||$ for $U_L$&' + '&'.join(str('%.3f' % avglgQgs[n]) for n in nrange) + '\\\\')
print('average $\\log_2||u^*||$ for $\\units\\LL$&' + '&'.join(str('%.3f' % avglgUgs[n]) for n in nrange) + '\\\\')
print('% average $||u^*||^2$ for $\\units\\LL$&' + '&'.join(str('%.3f' % avgsqUgs[n]) for n in nrange) + '\\\\')
print('average $\\log_2(\\#(\\units\\LL/U_L))$&' + '&'.join(str('%.3f' % avglgindex[n]) for n in nrange) + '\\\\')
