#!/bin/bash

root -q -b root://dtn-eic.jlab.org//work/eic2/EPIC/RECO/<current_file> -e 'for (auto b : *events->GetListOfLeaves()) { if (events->GetBranch(b->GetName()) == nullptr) continue; cout << events->GetBranch(b->GetName())->GetTotalSize() << " " << b->GetName() << endl; }' | sort -n > branch_size_current.txt
root -q -b root://dtn-eic.jlab.org//work/eic2/EPIC/RECO/<default_file> -e 'for (auto b : *events->GetListOfLeaves()) { if (events->GetBranch(b->GetName()) == nullptr) continue; cout << events->GetBranch(b->GetName())->GetTotalSize() << " " << b->GetName() << endl; }' | sort -n > branch_size_default.txt
python plot.py -c branch_size_current.txt -d branch_size_default.txt
