class StructureFileNotProvided(Exception):
    def __init__(self):
        self.message = "File provided is not StructureType(.st) file."
        super().__init__(self.message)
