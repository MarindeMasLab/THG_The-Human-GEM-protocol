import unittest

import pandas as pd
from functions import function_reac_identification as reactions


class ReactionIdentificationCharacterizationTests(unittest.TestCase):
    def test_jaccard_uses_list_lengths_for_union(self):
        self.assertEqual(reactions.jaccard(["a", "b"], ["b", "c"]), 1 / 3)
        self.assertEqual(reactions.jaccard(["a", "b"], ["a", "b"]), 1.0)

    def test_execute_jaccard_matches_same_direction_and_reverse_direction(self):
        left = pd.DataFrame(
            [
                ["left_same", None, None, None, None, "A,B", "C"],
                ["left_reverse", None, None, None, None, "X", "Y,Z"],
                ["left_nomatch", None, None, None, None, "Q", "R"],
            ]
        )
        right = pd.DataFrame(
            [
                [None, "right_same", None, None, None, "A,B", "C"],
                [None, "right_reverse", None, None, None, "Y,Z", "X"],
            ]
        )

        result = reactions.execute_jaccard(left, right)

        self.assertEqual(result[0].tolist(), ["left_same", "left_reverse"])
        self.assertEqual(result[3].tolist(), ["right_same", "right_reverse"])
        self.assertEqual(result[6].tolist(), [2.0, 2.0])

    def test_execute_jaccard_returns_empty_dataframe_with_expected_columns(self):
        left = pd.DataFrame([["left", None, None, None, None, "A", "B"]])
        right = pd.DataFrame([[None, "right", None, None, None, "C", "D"]])

        result = reactions.execute_jaccard(left, right)

        self.assertTrue(result.empty)
        self.assertEqual(list(result.columns), list(range(7)))

    def test_identify_reaction_extracts_current_space_delimited_format(self):
        reaction = """
        <reaction metaid="MAR0001" id="MAR0001" reversible="false">
          <annotation>
            <rdf:li rdf:resource="http://identifiers.org/rhea/R12345"/>
            <rdf:li rdf:resource="http://identifiers.org/ec-code/+1.2.3.4"/>
          </annotation>
          <listOfReactants>
            <speciesReference species="MAM02040c" stoichiometry="1"/>
            <speciesReference species="MAM00001c" stoichiometry="2"/>
          </listOfReactants>
          <listOfProducts>
            <speciesReference species="MAM02039c" stoichiometry="1"/>
            <speciesReference species="MAM00002m" stoichiometry="3"/>
          </listOfProducts>
          <geneProductAssociation>
            <geneProductRef geneProduct="GENE1"/>
          </geneProductAssociation>
        </reaction>
        """

        parsed = reactions.identify_reaction(
            reaction,
            7,
            r"MAM[0-9]+[a-z][0-9]*",
            "MAM02040",
            "MAM02039",
        )

        fields = parsed.strip().split(" ")

        self.assertEqual(fields[0:2], ["7", "R12345"])
        self.assertEqual(set(fields[2].split(",")), {"c", "m"})
        self.assertEqual(set(fields[3].split(",")), {"c", "m"})
        self.assertEqual(
            fields[4:],
            [
                "false",
                "MAM00001c",
                "MAM00002m",
                "2,3",
                "1,2,1,3",
                "1.2.3.4",
                "GENE1",
                "MAR0001",
            ],
        )


if __name__ == "__main__":
    unittest.main()
