# Validate mode shape
import pandas as pd
from acoustic_analyser import frame
import numpy as np

test_frame = frame.from_file(
    member_file="test_member.json", constraint_file="test_constraint.csv", debug=False
)

free_end = test_frame.fixed_end(member_id=0)
free_end = test_frame.free_end(member_id=3)

mode_shape, _, orignal_shape = test_frame.get_mode_shape(
    natural_freq=13.897, step_size=0.005, scaling_factor=0.25
)

ux = pd.read_csv("ux.csv", index_col=0)
uy = pd.read_csv("uy.csv", index_col=0)
ansys_result = pd.read_csv("nodes.csv", index_col=0)
ansys_result["ux_ansys"] = ux["ux"]
ansys_result["uy_ansys"] = uy["uy"]
ansys_result["y"] = ansys_result["y"].round(6)
ansys_result["x"] = ansys_result["x"].round(6)

ansys_result.drop(columns=["z"], inplace=True)

code_result = pd.DataFrame()
code_result["x"] = orignal_shape[:, 0].round(6)
code_result["y"] = orignal_shape[:, 1].round(6)
code_result["ux_python"] = mode_shape[:, 0] - orignal_shape[:, 0]
code_result["uy_python"] = mode_shape[:, 1] - orignal_shape[:, 1]

results = pd.merge(ansys_result, code_result, how="left", on=["x", "y"])

results.to_csv("mode_shape_comparision.csv")
