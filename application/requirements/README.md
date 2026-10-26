# Pinned environments

A file here named `<env>.txt` (`appenv.txt`, `funmapenv.txt`, `polyfunenv.txt`)
is installed with `pip install -r` and **takes precedence** over the loose
package list in `setup.sh`. That is what makes the Python environments
reproducible: without a file here, pip resolves to whatever is current, and
this analysis ran on `numpy 2.4.6 / pandas 3.0.6 / pyarrow 26.0.0 /
torch 2.14.1+cu130`, which a later install would not reproduce.

To create them, run `setup.sh` once and copy the frozen output:

    cp "$APP_ROOT"/pins/appenv.txt     application/requirements/appenv.txt
    cp "$APP_ROOT"/pins/funmapenv.txt  application/requirements/funmapenv.txt
    cp "$APP_ROOT"/pins/polyfunenv.txt application/requirements/polyfunenv.txt

then commit them. `setup.sh` writes `pins/<env>.txt` on every run, so the
difference between a pinned file and the current state is always visible with
`diff`.

Upstream repository commits live in `$APP_ROOT/pins/upstream.txt` and are
pinned by `FUNMAP_COMMIT` / `POLYFUN_COMMIT` / `FB_COMMIT` in `setup.sh`.
