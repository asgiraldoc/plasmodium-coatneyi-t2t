#! /usr/bin/env python3

import sys
import operator
import collections
import re

import utils.myFile
import utils.myTools
from utils.myTools import file
import utils.myGenomes
import utils.myPsOutput
import utils.myKaryoDrawer

__doc__ = """
    Compare two genomes using an orthologues list (ancGenes or 2 columns orthologues gene list)
        
        - can draw dotplot or karyotype
        - can print orthologous gene list, orthologous chromosome list, gene difference, orthologs count
        - can identify and count major rearrangements (fusions, translocations, inversions)
        
    
    Usage:    misc.compareGenomes.py genesST.Homo.sapiens.list.bz2 genesST.Mus.musculus.list.bz2 ancGenes.Euarchontoglires.list.bz2  > dotplot_Human_Mouse.ps
        misc.compareGenomes.py genesST.Homo.sapiens.list.bz2 genesST.Mus.musculus.list.bz2 ancGenes.Euarchontoglires.list.bz2 -mode=drawKaryotype -minChrSize=200 > Karyo_human_min200genes.ps
        misc.compareGenomes.py genesST.Homo.sapiens.list.bz2 genesST.Mus.musculus.list.bz2 ancGenes.Euarchontoglires.list.bz2 -mode=identifyRearrangements -rearr:minBlockSize=10
"""

# Palette of 23 distinguishable colors
COLORS = [
    "#1f77b4", "#aec7e8", "#ff7f0e", "#000000", "#2ca02c",
    "#98df8a", "#d62728", "#ff9896", "#9467bd", "#c5b0d5",
    "#8c564b", "#c49c94", "#e377c2", "#f7b6d2", "#7f7f7f",
    "#c7c7c7", "#bcbd22", "#dbdb8d", "#17becf", "#9edae5",
    "#ffbb78", "#555555", "#ffe118"
]

def hex_to_ps_color(hexcode):
    """
    Receives a string like "#RRGGBB" and returns a string "#R:G:B"
    with R, G, B as integers [0..255], allowing setColor to parse it.
    """
    r = int(hexcode[1:3], 16)
    g = int(hexcode[3:5], 16)
    b = int(hexcode[5:7], 16)
    return f"#{r}:{g}:{b}"

# ---------------------------------------------------
# Function for "natural" alphanumeric sorting
def natural_key(s):
    """
    Splits the string into chunks of digits and non-digits,
    converting digits to integers so that 10 > 2, etc.
    """
    return [int(text) if text.isdigit() else text.lower()
            for text in re.split(r'(\d+)', s) if text]
# ---------------------------------------------------

modes = ["drawMatrix", "drawKaryotype", "printOrthologuesList",
         "printOrthologuesCount", "printGeneDiff", "printOrthologousChrom",
         "identifyRearrangements"] # <-- ADDED NEW MODE

arguments = utils.myTools.checkArgs(
    [("studiedGenome", file), ("referenceGenome", file), ("orthologuesList", file)],
    [("includeGaps", bool, False), ("includeScaffolds", bool, True),
     ("includeRandoms",bool,False), ("includeNones",bool,False),
     ("reverse",bool,False),
     ("mode",str,modes),
     ("rearr:minBlockSize", int, 5), # <-- ADDED NEW PARAMETER
     ("orthoslist:fullgenenames",bool,False),
     ("orthoschr:minHomology",int,90),
     ("minChrSize",int,0),
     ("matrix:scaleY",bool,False), ("matrix:pointSize",float,-1),
     ("sortBySize",bool,False),
     ("matrix:colorFile",str,""), ("matrix:defaultColor",str,"black"),
     ("matrix:penColor",str,"black"),
     ("karyo:landscape",bool,False),
     ("ps:backgroundColor",str,"")],
    __doc__
)

# File loading
genesAnc = utils.myGenomes.Genome(arguments["orthologuesList"])
genome1 = utils.myGenomes.Genome(arguments["studiedGenome"], ancGenes=genesAnc)
genome2 = utils.myGenomes.Genome(arguments["referenceGenome"], ancGenes=genesAnc)
if arguments["reverse"]:
    genome1, genome2 = genome2, genome1

# Initial construction of contig lists
chr1 = list(genome1.chrList[utils.myGenomes.ContigType.Chromosome])
chr2 = list(genome2.chrList[utils.myGenomes.ContigType.Chromosome])
if arguments["includeScaffolds"]:
    chr1.extend(genome1.chrList[utils.myGenomes.ContigType.Scaffold])
    chr2.extend(genome2.chrList[utils.myGenomes.ContigType.Scaffold])
if arguments["includeRandoms"]:
    chr1.extend(genome1.chrList[utils.myGenomes.ContigType.Random])
    chr2.extend(genome2.chrList[utils.myGenomes.ContigType.Random])
if arguments["includeNones"]:
    chr1.extend(genome1.chrList[utils.myGenomes.ContigType.Unknown])
    chr2.extend(genome2.chrList[utils.myGenomes.ContigType.Unknown])

print(len(chr1), len(chr2), file=sys.stderr)

# Filter by minimum chromosome size
chr1 = [c for c in chr1 if len(genome1.lstGenes[c]) >= arguments["minChrSize"]]
chr2 = [c for c in chr2 if len(genome2.lstGenes[c]) >= arguments["minChrSize"]]

# Natural (alphanumeric) sorting of chromosomes
chr1.sort(key=natural_key)
chr2.sort(key=natural_key)

# Build ortholog tables
table12 = genome1.buildOrthosTable(chr1, genome2, chr2, arguments["includeGaps"], genesAnc)
table21 = genome2.buildOrthosTable(chr2, genome1, chr1, arguments["includeGaps"], genesAnc)

###########################################
def drawMatrix():
    print("Displaying ", end=' ', file=sys.stderr)

    if arguments["sortBySize"]:
        chr1.sort(key=lambda c: len(genome1.lstGenes[c]), reverse=True)
        chr2.sort(key=lambda c: len(genome2.lstGenes[c]), reverse=True)

    utils.myPsOutput.printPsHeader()
    if arguments["ps:backgroundColor"]:
        utils.myPsOutput.drawBox(0,0, 21,29.7,
                                 arguments["ps:backgroundColor"],
                                 arguments["ps:backgroundColor"])
    sys.stderr.write('.')

    colors_obj = None
    if arguments["matrix:colorFile"]:
        colors_obj = utils.myGenomes.Genome(arguments["matrix:colorFile"])

    # Scales
    nb = sum(len(table12[c]) for c in table12)
    scaleX = 19. / float(nb)
    scaleY = (19. / float(sum(len(table21[c]) for c in table21))
              if arguments["matrix:scaleY"] else scaleX)
    dp = scaleX if arguments["matrix:pointSize"] < 0 else arguments["matrix:pointSize"]
    sys.stderr.write('.')

    # Fixed color per chromosome of species 1
    color_chr1 = {c: COLORS[i % len(COLORS)] for i, c in enumerate(chr1)}

    def prepareGenome(dicOrthos, lst, func):
        i = 0
        y = 0
        lstNum = {}
        for c in lst:
            func(c, y, len(dicOrthos[c]))
            y += len(dicOrthos[c])
            for (gene,_) in dicOrthos[c]:
                lstNum[(c,gene)] = i
                i += 1
        func(None, y, None)
        return lstNum

    dl1 = float(sum(len(table21[c]) for c in table21)) * scaleY
    def line1(c, x, l):
        utils.myPsOutput.drawLine(1 + x*scaleX, 1, 0, dl1,
                                  arguments["matrix:penColor"])
        if c:
            utils.myPsOutput.drawText(1 + (x+l/2)*scaleX, 0.7, c,
                                      arguments["matrix:penColor"])

    def line2(c, x, l):
        utils.myPsOutput.drawLine(1, 1 + x*scaleY, 19, 0,
                                  arguments["matrix:penColor"])
        if c:
            print("90 rotate")
            utils.myPsOutput.drawText(1 + (x+l/2)*scaleY, -0.9, c,
                                      arguments["matrix:penColor"])
            print("-90 rotate")

    lstNum1 = prepareGenome(table12, chr1, line1)
    sys.stderr.write('.')
    lstNum2 = prepareGenome(table21, chr2, line2)
    sys.stderr.write('.')

    print("0 setlinewidth")
    for c1 in table12:
        for (i1,t) in table12[c1]:
            xx = 1 + float(lstNum1[(c1,i1)]) * scaleX
            for (c2,i2) in t:
                coul = color_chr1[c1]
                if colors_obj is not None:
                    tmp = set(colors_obj.getPosition(
                              genome1.lstGenes[c1][i1].names
                              + genome2.lstGenes[c2][i2].names))
                    for (cc,ii) in genesAnc.getPosition(
                              genome1.lstGenes[c1][i1].names
                              + genome2.lstGenes[c2][i2].names):
                        tmp.update(colors_obj.getPosition(
                                   genesAnc.lstGenes[cc][ii].names))
                    if tmp:
                        coul = tmp.pop()[0]
                yy = 1 + lstNum2[(c2,i2)]*scaleY
                utils.myPsOutput.drawBox(xx, yy, dp, dp, coul, coul)

    # Labels
    utils.myPsOutput.drawText(4, 0.3,
                              arguments["referenceGenome"]
                              if arguments["reverse"]
                              else arguments["studiedGenome"],
                              arguments["matrix:penColor"])
    print("90 rotate")
    utils.myPsOutput.drawText(4, -0.5,
                              arguments["studiedGenome"]
                              if arguments["reverse"]
                              else arguments["referenceGenome"],
                              arguments["matrix:penColor"])
    print("-90 rotate")
    utils.myPsOutput.printPsFooter()
    print(" OK", file=sys.stderr)


##############################################################################################
def drawKaryotype():
    # Header and background
    (lx,ly) = utils.myPsOutput.printPsHeader(arguments["karyo:landscape"])
    if arguments["ps:backgroundColor"]:
        utils.myPsOutput.drawBox(0,0, lx,ly,
                                 arguments["ps:backgroundColor"],
                                 arguments["ps:backgroundColor"])

    # Prepare data: list of tuples (chromosome, [orthologs...])
    data = []
    for c in chr1:
        newl = []
        for (_,val) in table12.get(c, []):
            newl.append(val[0][0] if val else None)
        data.append((c, newl))

    # Color mapping: using the natural order already applied to chr2
    color_chrom = {}
    for i, c in enumerate(chr2):
        hexcode = COLORS[i % len(COLORS)]
        color_chrom[c] = hex_to_ps_color(hexcode)

    print("Displaying ...", end=' ', file=sys.stderr)
    utils.myKaryoDrawer.drawKaryo(
        data,
        arguments,
        x0=1, y0=1,
        lx=lx-2,
        ly=ly-2,
        bysize=arguments["sortBySize"],
        color_chrom=color_chrom
    )
    utils.myPsOutput.printPsFooter()
    print("OK", file=sys.stderr)


#########################################################################
def printOrthologuesCount():
    print(utils.myFile.myTSV.printLine([""] + chr2))
    for c1 in chr1:
        count = collections.defaultdict(int)
        for (_,t) in table12[c1]:
            for (c2,_) in t:
                count[c2] += 1
        print(utils.myFile.myTSV.printLine([c1] + [count[c2] for c2 in chr2]))


###############################################################################
def printOrthologuesList():
    def printGene(g):
        s = list(g)
        s[-1] = "/".join(s[-1]) if arguments["orthoslist:fullgenenames"] else s[-1][0]
        return s

    for c1 in chr1:
        for (i1,t) in sorted(table12[c1]):
            g1 = genome1.lstGenes[c1][i1]
            for (c2,i2) in sorted(t):
                print(utils.myFile.myTSV.printLine(
                    printGene(g1) + printGene(genome2.lstGenes[c2][i2])
                ))


####################################################
def printGeneDiff():
    def getGeneTxt(g):
        return "/".join(g.names) + ":%s:%d-%d:%d" % g[:4]

    all_links = set()
    combin = utils.myTools.myCombinator()
    for c1 in table12:
        for (i1,t) in table12[c1]:
            combin.addLink([(1,c1,i1)] + [(2,c2,i2) for (c2,i2) in t])
    for c2 in table21:
        for (i2,t) in table21[c2]:
            combin.addLink([(2,c2,i2)] + [(1,c1,i1) for (c1,i1) in t])
    for g in combin:
        e1 = [getGeneTxt(genome1.lstGenes[c][i]) for (x,c,i) in g if x == 1]
        e2 = [getGeneTxt(genome2.lstGenes[c][i]) for (x,c,i) in g if x == 2]
        if not e1:
            print("+", end=' ')
        elif not e2:
            print("-", end=' ')
        elif len(e1)==1 and len(e2)==1:
            print("=", end=' ')
        elif len(e1)>1 and len(e2)==1:
            print("--", end=' ')
        elif len(e1)==1 and len(e2)>1:
            print("++", end=' ')
        else:
            print("**", end=' ')
        print(" ".join(e1 + e2))


#############################################################
def printOrthologousChrom():
    for c1 in chr1:
        count = collections.defaultdict(int)
        for (_,t) in table12[c1]:
            for (c2,_) in t:
                count[c2] += 1
        res = [c1]
        t = sorted(count.items(), key=operator.itemgetter(1))
        n = (sum(count.values()) * arguments["orthoschr:minHomology"]) / 100
        while n > 0 and t:
            x = t.pop()
            res.append(f"{x[0]} ({x[1]})")
            n -= x[1]
        print(utils.myFile.myTSV.printLine(res))

#########################################################################
# NEW FUNCTION ADDED
#########################################################################
def identifyRearrangements():
    """
    Identifies and counts fusions, translocations, and inversions
    based on the collinearity of orthologous gene blocks.
    """
    print("[INFO] Starting structural rearrangements analysis...", file=sys.stderr)
    min_block_size = arguments["rearr:minBlockSize"]
    
    # 1. Build synteny blocks
    # Block format: (chr2, strand, num_genes, start_g1, end_g1, start_g2, end_g2)
    synteny_blocks = collections.defaultdict(list)

    for c1 in chr1:
        if not table12.get(c1):
            continue

        # Sort genes by their position in genome 1
        sorted_genes_g1 = sorted(table12[c1], key=lambda x: x[0])
        
        current_block = []
        last_c2 = None
        last_g2_idx = None
        # Orientation (strand) is initially determined with the first pair of genes
        last_strand = 0

        for g1_idx, orthologs in sorted_genes_g1:
            if not orthologs:
                continue
            
            # Use the first ortholog as representative for simplicity
            c2, g2_idx = orthologs[0]

            # Determine strand orientation change
            strand = 0
            if last_g2_idx is not None:
                if g2_idx > last_g2_idx:
                    strand = 1 # Positive strand (+)
                elif g2_idx < last_g2_idx:
                    strand = -1 # Negative strand (-) (possible inversion)
            
            # A block ends if the destination chromosome (c2) or orientation changes
            if c2 != last_c2 or (strand != 0 and last_strand != 0 and strand != last_strand):
                if len(current_block) >= min_block_size:
                    # Save the complete block
                    start_g1 = current_block[0][0]
                    end_g1 = current_block[-1][0]
                    g2_indices = sorted([b[1] for b in current_block])
                    synteny_blocks[c1].append( (last_c2, last_strand, len(current_block), start_g1, end_g1, g2_indices[0], g2_indices[-1]) )
                
                # Start a new block
                current_block = []
                last_strand = strand
            
            current_block.append((g1_idx, g2_idx))
            last_c2 = c2
            last_g2_idx = g2_idx
            if last_strand == 0: # Assign orientation if it is the start of a series
                last_strand = strand

        # Save the last block in the list if it is large enough
        if len(current_block) >= min_block_size:
            start_g1 = current_block[0][0]
            end_g1 = current_block[-1][0]
            g2_indices = sorted([b[1] for b in current_block])
            synteny_blocks[c1].append( (last_c2, last_strand, len(current_block), start_g1, end_g1, g2_indices[0], g2_indices[-1]) )

    print(f"[INFO] Identified {sum(len(b) for b in synteny_blocks.values())} syntenic blocks (min. {min_block_size} genes).", file=sys.stderr)

    # 2. Analyze the blocks to count events
    fusions = 0
    translocations = 0
    inversions = 0
    
    print("\n--- Structural Rearrangements Report ---")
    
    for c1, blocks in synteny_blocks.items():
        if not blocks:
            continue
            
        # Count source chromosomes (c2) for this c1
        source_chroms_count = collections.Counter(b[0] for b in blocks)
        
        # --- Fusion/Fission Detection ---
        # If a chromosome has significant blocks from more than one source chromosome
        if len(source_chroms_count) > 1:
            fusions += 1
            print(f"\n[FUSION/FISSION DETECTED] Chromosome '{c1}' in {genome1.name} contains blocks from {len(source_chroms_count)} chromosomes of {genome2.name}:")
            for c2, count in source_chroms_count.items():
                print(f"  - Source: '{c2}' ({count} block{'s' if count > 1 else ''})")
        
        # If there is more than one block, we can analyze translocations and inversions
        if len(blocks) > 1:
            # Sort blocks by their start position in genome 1
            sorted_blocks_on_c1 = sorted(blocks, key=lambda x: x[3])
            major_source_chrom = source_chroms_count.most_common(1)[0][0]
            
            for i, block in enumerate(sorted_blocks_on_c1):
                c2, strand, num_genes, start_g1, _, _, _ = block
                
                # --- Interchromosomal Translocation Detection ---
                if c2 != major_source_chrom:
                    translocations += 1
                    print(f"\n[TRANSLOCATION] Block from '{c2}' ({num_genes} genes) found in '{c1}' (whose major source chromosome is '{major_source_chrom}').")

                # --- Inversion Detection ---
                # Considered an inversion if the block is on the negative strand
                if strand == -1:
                    inversions += 1
                    print(f"\n[INVERSION] Block in '{c1}' (of {num_genes} genes) mapping to '{c2}' is in inverted orientation.")
    
    print("\n-----------------------\n--- Final Summary ---\n-----------------------")
    print(f"Total Fusion/Fission events detected: {fusions}")
    print(f"Total Interchromosomal Translocations detected: {translocations}")
    print(f"Total Inversions detected: {inversions}")
    print("\nNote: Event counts are an approximation based on the grouping of syntenic blocks.")


# Execute the requested mode
eval(str(arguments["mode"]))()
