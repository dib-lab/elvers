from io import StringIO
import sys


class CaptureStdout(list):
    """
    A utility object that uses a context manager
    to capture stdout from Snakemake. Useful when
    creating the directed acyclic graph.
    """
    def __init__(self,passthru=False):
        # Boolean: should we pass everything through to stdout?
        # (this object is only functional if passthru is False)
        self.passthru = passthru

    def __enter__(self):
        # Open a new context with this CaptureStdout
        # object. This happens when we say
        # "with CatpureStdout('label') as output:"

        # If we are just passing input on to output, pass thru
        if self.passthru:
            return self

        # Otherwise, we want to swap out sys.stdout with
        # a StringIO object that will save stdout.
        # 
        # Save the existing stdout object so we can
        # restore it when we're done
        self._stdout = sys.stdout
        # Now swap out stdout 
        sys.stdout = self._stringio = StringIO()
        return self

    def __exit__(self, *args):
        # If we are just passing input on to output, pass thru
        if self.passthru:
            return self

        # This entire class extends the list class,
        # so we call self.extend() to add a list to 
        # the end of self (in this case, all the new
        # lines from our StringIO object).
        self.extend(self._stringio.getvalue().splitlines())
        del self._stringio    # free up some memory

        # Clean up by setting sys.stdout back to what
        # it was before we opened up this context.
        sys.stdout = self._stdout

