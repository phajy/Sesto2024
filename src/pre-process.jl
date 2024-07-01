# assume user has `xz` and `tar` installed

REFLIONX_DATA_DIR = "./data/reflionx"
run(`mkdir -p $REFLIONX_DATA_DIR`)
run(`tar xf ./data/grid.tar.xz -C $REFLIONX_DATA_DIR`)
