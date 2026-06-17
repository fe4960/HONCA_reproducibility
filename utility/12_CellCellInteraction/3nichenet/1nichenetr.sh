#!/bin/sh

source ~/.condainit

conda activate archr

R < HCA_ON/scripts/12_CellCellInteraction/3nichenet/1nichenetr.R --no-save > HCA_ON/scripts/12_CellCellInteraction/3nichenet/1nichenetr.out 2> HCA_ON/scripts/12_CellCellInteraction/3nichenet/1nichenetr.err
