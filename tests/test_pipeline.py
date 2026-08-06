import pandas as pd

from staging.my_source.pipeline import normalize_negative_control_values


def test_negative_control_value_moves_to_moa():
    bio_df = pd.DataFrame(
        {
            "moa": ["", "Inhibitor"],
            "value": ["Negative Control", 12.5],
        }
    )

    result = normalize_negative_control_values(bio_df)

    assert result.loc[0, "moa"] == "Negative control"
    assert pd.isna(result.loc[0, "value"])
    assert result.loc[1, "moa"] == "Inhibitor"
    assert result.loc[1, "value"] == 12.5
