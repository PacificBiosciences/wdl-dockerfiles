"""
Patch truvari.phab's reference-sequence extraction to write directly via
`samtools faidx -o <file>` instead of relying on pysam's captured-stdout
return value.

pysam 0.24.0 (the only pysam release with Python 3.14 support -- it's the
first release published after Python 3.14 existed, so every earlier
release predates and cannot build against it) silently stopped capturing
`samtools faidx`'s stdout: pysam's dispatch wrapper only calls
`samtools_set_stdout_fn()` for methods listed in its own MAP_STDOUT_OPTIONS
dict (view/mpileup/depad/calmd), and faidx isn't one of them. The command's
real output still goes to the process's actual stdout, but the string
pysam returns to the caller is empty, so truvari.phab.VCFtoHaplotypes.set_regions
writes an empty file and the immediately following `samtools faidx` on that
empty file fails with "File truncated at line 1".
"""
import truvari.phab as p

path = p.__file__
with open(path) as f:
    src = f.read()

old = '''        # Pull sequences
        out_fn = truvari.make_temp_filename(suffix='.fa')
        with open(out_fn, 'w') as fout:
            fout.write(samtools.faidx(
                self.reference_fn, "-r", regions_file_name))
        # Facilitate fetching
        samtools.faidx(out_fn)
        self.ref_haps_fn = out_fn'''

new = '''        # Pull sequences
        out_fn = truvari.make_temp_filename(suffix='.fa')
        samtools.faidx(self.reference_fn, "-r", regions_file_name, "-o", out_fn)
        # Facilitate fetching
        samtools.faidx(out_fn)
        self.ref_haps_fn = out_fn'''

assert old in src, (
    f"expected phab.py block not found in {path} -- "
    "truvari source may have changed since this patch was written"
)
src = src.replace(old, new, 1)
with open(path, 'w') as f:
    f.write(src)
print(f"Patched {path}")
