
class WeakMatroid:
    def __init__(self, matroid) -> None:
        self.matroid = matroid

    def __eq__(self, __value: object) -> bool:
        """
        In our work, we can work with "uniqueness up to
        isomorphism". Hence, two matroids are equal here
        if they are isomorphic to one another.
        """
        return self.matroid.is_isomorphic(__value.matroid)
    
    def __ne__(self, __value: object) -> bool:
        return not self.__eq__(__value)
    
    def __hash__(self) -> int:
        # Sage hashes basis matroids the following way:
        # hash((matroid.groundset(), matroid.bases_count(), matroid._weak_invariant()))
        # And since this is the type of matroid we deal with here,
        # it's sufficient to use like so.
        return hash(self.matroid)
    
    def __str__(self) -> str:
        return str(self.matroid)

