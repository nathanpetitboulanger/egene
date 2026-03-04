sequence = "AGTCGTAC"

revert_dict = {
    "A": "T",
    "T": "A",
    "C": "G",
    "G": "C",
}


class DNA_manipulator:
    def __init__(self, sequence):
        self.sequence = sequence
        self.comp_dict = {
            "A": "T",
            "T": "A",
            "C": "G",
            "G": "C",
        }

    @property
    def list_seq(self) -> list[str]:
        return [letter for letter in self.sequence]

    @property
    def revert_list(self) -> list[str]:
        return self.list_seq[::-1]

    @property
    def comp_list(self) -> list[str]:
        return [self.comp_dict[key] for key in self.list_seq]

    @property
    def comp_seq(self) -> str:
        return "".join(self.comp_list)

    @property
    def revers_seq(self) -> str:
        return "".join(self.revert_list)

    @property
    def revers_comp_list(self) -> list[str]:
        return [self.comp_dict[key] for key in self.comp_list]

    @property
    def revers_comp_seq(self) -> str:
        return "".join(self.revers_comp_list)

    pass


seq = DNA_manipulator(sequence)

seq.revers_comp_seq
