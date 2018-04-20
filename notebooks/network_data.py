
"""
Utility data manipulation routines for use with the network data.
"""

import numpy as np

tag_filename = "tagged_network.tsv"
network_file = "reg_net.csv"
heatmap_filename = "heatmap.txt"
fixed_heatmap_filename = "tab_first_heatmap.txt"
default_colname = "pairwise_split_499"
coding_filename = "coding.txt"

class GeneCategories(object):

    def __init__(self, tag_filename=tag_filename):
        self.tag_filename = tag_filename
        self.read_data()

    def read_data(self):
        f = open(self.tag_filename)
        L = []
        while 1:
            ln = f.readline()
            if not ln:
                break
            L.append(ln[:-1].split("\t"))
        A = self.matrix = np.array(L, dtype=object)
        id_to_name = {}
        for row in A:
            id_to_name[row[0]] = row[2]
        self.id_to_name = id_to_name
        self.identities = A[:, 0]
        self.names = A[:, 2]

    def fix_heatmap(self, heatmap_filename=heatmap_filename, fixed_heatmap_filename=fixed_heatmap_filename):
        "Change id's to gene names in heatmap file and fix other formatting issues."
        
        f = open(coding_filename)
        sample_map = {}
        headers = f.readline().strip().split("\t")
        # actually
        headers = "TIME DONOR TREATMENT CODE SAMPLE".split()
        #print "headers " + repr(headers)
        while 1:
            line = f.readline().strip()
            if not line:
                break
            mapping = line.split("\t")
            #print "mapping " + repr(mapping)
            D = dict(zip(headers, mapping))
            sample_map[D["SAMPLE"]] = D
            #print D
            #break
        def condition_description(sample_id):
            D = sample_map.get(sample_id)
            if D is None:
                return sample_id
            return ".".join([D["DONOR"], D["TIME"], D["TREATMENT"], D["CODE"] ])
        f = open(heatmap_filename)
        outf = open(fixed_heatmap_filename, "w")
        #outf.write(f.read())
        headers = f.readline()
        if headers[0] != "\t":
            outf.write("\t")  # add a leading tab
        headernames = headers.split()
        descriptors = map(condition_description, headernames)
        outf.write("\t".join(descriptors) + "\n")
        #outf.write(headers)
        while 1:
            line = f.readline()
            if not line:
                break
            (identifier, etcetera) = line.split("\t", 1)
            outf.write(self.mapid(identifier))
            outf.write("\t")
            outf.write(etcetera)
        outf.close()

    def mapid(self, identity):
        return self.id_to_name.get(identity, identity)

    def map_csv_gene_network(self, network_file=network_file):
        f = open(network_file)
        def get_entries():
            txt = f.readline().strip()
            return txt.split(",")
        headers = get_entries()
        assert len(headers) >= 2, "expect source and target columns in csv file " + repr(from_file)
        result = []
        while 1:
            entries = get_entries()
            if len(entries) == 1 and len(entries[0]) == 0:
                break
            source = entries[0]
            target = entries[1]
            result.append((self.mapid(source), self.mapid(target)))
        return result

    def network_widget(self, verbose=True, N=None):
        from jp_gene_viz import dNetwork, getData, dGraph, grid_forest
        dNetwork.load_javascript_support()
        if verbose:
            print("Reading network graph")
        G = dGraph.WGraph()
        edges = self.map_csv_gene_network()
        for (regulator, target) in edges:
            G.add_edge(regulator, target, 1.0)
        if verbose:
            print("Assigning grid layout")
        (layout, rectangles) = grid_forest.forest_layout(G)
        if verbose:
            print("Initializing network widget")
        if N is None:
            N = dNetwork.NetworkDisplay()
        self.network = N
        N.load_data(G, layout, draw=False)
        if verbose:
            print("Configuring network widget parameters.")
        N.container_dropdown.value = dNetwork.CANVAS
        #N.node_categories = self.top_categorization()
        self.apply_categorization(default_colname, N)
        self.apply_colorization(default_colname, N)
        N.label_rectangles = True
        N.rectangle_color = "#999999"
        return N

    def colorize_shortest_path(self, from_gene, to_gene, color="red", directed=True):
        #print "connecting", from_gene, to_gene
        N = self.network
        dG = N.display_graph
        edges = dG.edge_weights
        nodes = dG.node_weights
        back_path = {}
        adjacent = {node: set() for node in nodes}
        for (tf, gene) in edges:
            adjacent[tf].add(gene)
            if not directed:
                adjacent[gene].add(tf)
        horizon = {from_gene}
        # Point back from horizon nodes to previous node on shortest paths.
        # If there are many paths of the same length the algorithm choses one arbitrarily.
        while horizon and (to_gene not in back_path):
            #print "horizon", horizon
            #print "back", back_path
            next_horizon = set()
            for choose_gene in horizon:
                for next_gene in adjacent[choose_gene]:
                    if next_gene not in back_path and next_gene!=from_gene:
                        back_path[next_gene] = choose_gene
                        next_horizon.add(next_gene)
                del adjacent[choose_gene]  # prevent infinite loop if buggy
            horizon = next_horizon
        if to_gene not in back_path:
            raise KeyError("no path found from %s to %s" % (repr(from_tf), repr(to_gene)))
        gene = to_gene
        path = set()
        # walk backwards to construct the shortest path
        while gene != from_gene:
            back_gene = back_path[gene]
            del back_path[gene]  # prevents infinite loop if buggy
            path.add((back_gene, gene))
            gene = back_gene
        #print "path", from_gene, to_gene, path
        # path is now the shortest path: colorize it
        for (gene1, gene2) in path:
            if (gene1, gene2) in edges:
                svg_name = dG.edge_name(gene1, gene2)
            else:
                assert (gene2, gene1) in edges, "bad edge in path " + repr((gene1, gene2))
                svg_name = dG.edge_name(gene2, gene1)
            N.color_overrides[svg_name] = color
        N.draw()

    def apply_categorization(self, colname_pattern=default_colname, N=None):
        colname = self.match_colname(colname_pattern)
        if N is None:
            N = self.network
        N.node_categories = self.top_categorization(colname)

    def apply_colorization(self, colname_pattern=default_colname, N=None):
        colname = self.match_colname(colname_pattern)
        if N is None:
            N = self.network
        colorize = self.colorization(name=colname)
        for name in colorize:
            color = colorize[name]
            # Hack: refactor this!
            svg_name = N.display_graph.node_name(name)
            N.color_overrides[svg_name] = color
    
    def match_colname(self, colname_pattern):
        import fnmatch
        headers = self.matrix[0]
        matches = fnmatch.filter(headers, colname_pattern)
        nmatches = len(matches)
        if nmatches < 1:
            raise KeyError("No match found for %s" % repr(colname_pattern))
        if nmatches > 1:
            raise KeyError("Too many matches for %s: %s" % (repr(colname_pattern), matches))
        return matches[0]

    def network_and_heatmap_widget(self, verbose=True):
        from jp_gene_viz import LExpression
        L = LExpression.LinkedExpressionNetwork()
        self.network_and_heatmap = L
        self.network_widget(verbose=True, N=L.network)
        self.fix_heatmap()
        L.load_heatmap(fixed_heatmap_filename)
        return L
        
    def get_column_by_name(self, name=default_colname):
        result = None
        A = self.matrix
        (_, cols) = A.shape
        for c in range(cols):
            column = A[:, c]
            if column[0] == name:
                result = column
        if result is None:
            raise KeyError("column not found " + repr(name))
        return result

    def categories(self, name=default_colname):
        column = self.get_column_by_name(name)
        names = self.names
        result = {}
        for index in range(1, len(names)):
            result[names[index]] = column[index]
        return result
    
    def category_counts(self, name=default_colname):
        cats = self.categories(name)
        result = {}
        for name in cats:
            cat = cats[name]
            result[cat] = result.get(cat, 0) + 1
        return result

    def top_categories(self, name=default_colname, n=10):
        counts = self.category_counts(name)
        count_item = [(counts[cat], cat) for cat in counts]
        scounts = list(reversed(sorted(count_item)))
        return [x[1] for x in scounts[:n]]

    # first color is default
    colors = color_list_499 = ["lightgray", '#DA8DF0', '#1AB51A', 
        '#02CCDE', '#062833',  '#FA8E00', '#FFF527', '#5EAB9F', '#FA38A1', '#DB1002', '#0700DB']

    def colorization(self, name=default_colname):
        top_cats = self.top_categories(name)
        print "top categories are", top_cats
        cat_to_color = {}
        default_color = self.colors[0]
        for (index, cat) in enumerate(top_cats):
            cat_to_color[cat] = self.colors[index + 1]  # skip the default
        cats = self.categories(name)
        result = {}
        for n in cats:
            cat = cats[n]
            result[n] = cat_to_color.get(cat, default_color)
        return result

    def top_categorization(self, name=default_colname, by_ordinal=True):
        top_cats = self.top_categories(name)
        cat_to_ordinal = {cat: index for (index, cat) in enumerate(top_cats)}
        name_to_top_cat = {}
        cats = self.categories(name)
        for n in cats:
            cat = cats[n]
            if cat in top_cats:
                if by_ordinal:
                    name_to_top_cat[n] = cat_to_ordinal[cat]
                else:
                    name_to_top_cat[n] = cat
        return name_to_top_cat
