import ROOT
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("filename", type=str)
parser.add_argument("replace_old", type=str)
parser.add_argument("replace_new", type=str)
parser.add_argument("--skipSyst", default=False, action='store_true')

args = parser.parse_args()

f = ROOT.TFile(args.filename)
newf = ROOT.TFile(args.filename.replace(".root", "_corr.root"), "recreate")

d = f.GetDirectory("tagsDumper/trees")

for i, k in enumerate(d.GetListOfKeys()):
  #if "ggh_13TeV_UntaggedTag_0" not in k.GetName(): continue #DEBUG
  print(k)

  if args.skipSyst:
    if "sigma" in k.GetName(): continue

  t = f.Get(f"tagsDumper/trees/{k.GetName()}");
  old_name = t.GetName()
  t.SetName(t.GetName().replace(args.replace_old, args.replace_new))

  newf.cd()
  if i == 0:
    newf.mkdir("tagsDumper")
    newf.mkdir("tagsDumper/trees")
  newf.cd("tagsDumper/trees")
  newd = newf.GetDirectory("tagsDumper/trees")

  newt = t.CloneTree(-1, "fast");
  newt.SetDirectory(newd)

  newf.Write();

newf.Close()
