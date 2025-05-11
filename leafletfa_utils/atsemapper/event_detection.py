#!/usr/bin/env python3
import random
import logging
import networkx as nx # type: ignore
from typing import Dict, List, Tuple
from collections import defaultdict
from tqdm import tqdm # type: ignore
import gzip
import pandas as pd

class ATSEAnalyzer:
    def __init__(self):
        self.events = {}
        
    def build_splice_graph(self, junctions: Dict[str, Dict]):
        gene_graphs = {}
        gene_groups = defaultdict(dict)

        # Track statistics
        stats = {
            'total_junctions': len(junctions),
            'included_junctions': 0,
            'excluded_no_gene': 0,
            'junctions_per_gene': defaultdict(int),
            'single_junction_genes': set()  # Track genes with only one junction
        }

        # First pass - group junctions by gene
        for j_id, j_data in junctions.items():
            if not j_data['gene_ids']:
                stats['excluded_no_gene'] += 1
                continue
            gene_id = j_data['gene_ids'][0] # # Only takes the first gene ID
            gene_groups[gene_id][j_id] = j_data
            stats['included_junctions'] += 1

        # Build graphs and track single-junction genes
        for gene_id, gene_juncs in gene_groups.items():
            G = nx.Graph()

            for j_id, j_data in gene_juncs.items():
                if j_data['strand'] == '+':
                    donor = (j_data['chrom'], j_data['start'], 'donor')
                    acceptor = (j_data['chrom'], j_data['end'], 'acceptor')
                else:
                    donor = (j_data['chrom'], j_data['end'], 'donor')
                    acceptor = (j_data['chrom'], j_data['start'], 'acceptor')

                G.add_edge(donor, acceptor,
                          junction_id=j_id,
                          strand=j_data['strand'],
                          score=j_data['total_score'])

            gene_graphs[gene_id] = G
            num_junctions = len(gene_juncs)
            stats['junctions_per_gene'][gene_id] = num_junctions

            if num_junctions == 1:
                stats['single_junction_genes'].add(gene_id)

        # Calculate additional statistics
        stats['num_genes'] = len(gene_graphs)
        stats['num_single_junction_genes'] = len(stats['single_junction_genes'])
        stats['num_analyzable_genes'] = stats['num_genes'] - stats['num_single_junction_genes']
        stats['single_junction_pairs'] = stats['num_single_junction_genes']
        stats['avg_junctions_per_gene'] = stats['included_junctions'] / len(gene_graphs) if gene_graphs else 0
        stats['max_junctions_in_gene'] = max(stats['junctions_per_gene'].values()) if stats['junctions_per_gene'] else 0
        stats['min_junctions_in_gene'] = min(stats['junctions_per_gene'].values()) if stats['junctions_per_gene'] else 0

        # Print summary automatically
        print(f"""Splice Graph Building Summary:
    Total junctions: {stats['total_junctions']}
    - Included in graphs: {stats['included_junctions']}
    - Excluded (no gene): {stats['excluded_no_gene']}

    Gene Statistics:
    - Total genes: {stats['num_genes']}
    - Genes with single junction (can't analyze): {stats['num_single_junction_genes']} ({stats['num_single_junction_genes']/stats['num_genes']*100:.1f}% of genes)
    - Analyzable genes (>1 junction): {stats['num_analyzable_genes']} ({stats['num_analyzable_genes']/stats['num_genes']*100:.1f}% of genes)

    Junction Distribution:
    - Average junctions per gene: {stats['avg_junctions_per_gene']:.1f}
    - Max junctions in a gene: {stats['max_junctions_in_gene']}
    - Min junctions in a gene: {stats['min_junctions_in_gene']}""")

        return gene_graphs, stats

    def find_connected_junctions(self, G: nx.Graph, start_junction_id: str, visited: set) -> set:
        """
        Find all junctions connected by shared splice sites, including indirect connections.
        Returns empty set if only one junction is found.
        """
        # Get start junction's nodes
        start_nodes = []
        for u, v, data in G.edges(data=True):
            if data['junction_id'] == start_junction_id:
                start_nodes = [u, v]
                break
            
        if not start_nodes:
            return set()

        # Find all nodes in the connected component
        component_nodes = set(nx.node_connected_component(G, start_nodes[0]))

        # Get all junctions in this component
        connected = set()
        for u, v, data in G.edges(data=True):
            if u in component_nodes and v in component_nodes:
                junction_id = data['junction_id']
                if junction_id not in visited:
                    connected.add(junction_id)

        return connected if len(connected) > 1 else set()

    def analyze_splice_sites(self, G: nx.Graph) -> Dict:
        """
        Analyze splice site connectivity in a graph.
        """
        stats = {
            'num_nodes': G.number_of_nodes(),
            'num_edges': G.number_of_edges(),
            'donor_sites': len([n for n in G.nodes() if n[2] == 'donor']),
            'acceptor_sites': len([n for n in G.nodes() if n[2] == 'acceptor'])
        }

        # Analyze node degrees (how many junctions share each splice site)
        degrees = [d for _, d in G.degree()]
        stats['max_degree'] = max(degrees) if degrees else 0
        stats['isolated_sites'] = degrees.count(1)

        return stats

    def analyze_singleton_junctions(self, G: nx.Graph, junction_id: str) -> Dict:
        """Analyze a singleton junction for its properties."""
        junction_info = {}

        for u, v, data in G.edges(data=True):
            if data['junction_id'] == junction_id:
                junction_info = {
                    'junction_id': junction_id,
                    'donor_site': u[1],
                    'acceptor_site': v[1],
                    'strand': data['strand'],
                    'score': data['score'],  # number of reads
                    'chromosome': u[0]
                }
                break

        return junction_info
        
    def find_atse_groups(self, graphs: Dict[str, nx.Graph], min_splice_site_usage: float = 0.01) -> Tuple[Dict[str, Dict], Dict]:
        """
        Find alternative transcript splicing events (ATSEs) with splice site usage filtering.
        Returns dictionary with renamed ATSEs (gene_id_atse_1, etc.) and event counts.
        """
        # Use temporary naming during discovery, will rename later
        temp_atse_groups = {}
        event_counter = 0
        filtered_junctions = 0  # Track filtered junctions

        # Statistics tracking
        stats = {
            'total_junctions': sum(G.number_of_edges() for G in graphs.values()),
            'analyzable_junctions': 0,
            'junctions_in_atses': set(),
            'junction_counts': defaultdict(int),
            'singleton_junctions': [],  # Store info about singleton junctions
            'filtered_junctions': [],   # Track junctions filtered due to low splice site usage
            'split_atses': 0,          # Track how many original ATSEs were split
            'component_breakdown': defaultdict(int)  # Track component sizes after splitting
        }

        # First count analyzable junctions
        for gene_id, G in graphs.items():
            num_junctions = G.number_of_edges()
            if num_junctions >= 2:
                stats['analyzable_junctions'] += num_junctions

        # Find ATSEs
        for gene_id, G in graphs.items():
            visited = set()

            for _, _, data in G.edges(data=True):
                junction_id = data['junction_id']
                if junction_id not in visited:
                    connected = self.find_connected_junctions(G, junction_id, visited)

                    if connected:
                        # Get all splice sites in this connected component
                        splice_sites = set()
                        junction_data = {}

                        # First pass: collect all splice sites and junction data
                        for j_id in connected:
                            for u, v, data in G.edges(data=True):
                                if data['junction_id'] == j_id:
                                    splice_sites.add(u)
                                    splice_sites.add(v)
                                    junction_data[j_id] = {
                                    'donor': u,
                                    'acceptor': v,
                                    'score': data['score'],
                                    'strand': data['strand']
                                    }

                        # Calculate total reads at each splice site
                        site_total_reads = defaultdict(int)
                        for j_id, j_data in junction_data.items():
                            site_total_reads[j_data['donor']] += j_data['score']
                            site_total_reads[j_data['acceptor']] += j_data['score']

                        # Calculate usage proportions and filter junctions with low splice site usage
                        filtered_out = set()
                        junction_usage = {}

                        for j_id, j_data in junction_data.items():
                            donor_usage = j_data['score'] / site_total_reads[j_data['donor']]
                            acceptor_usage = j_data['score'] / site_total_reads[j_data['acceptor']]

                            # store the usage values for each junction 
                            if j_data['strand'] == '+':
                                # For positive strand, donor is 5' and acceptor is 3'
                                five_prime_usage = donor_usage
                                three_prime_usage = acceptor_usage
                            else:
                                # For negative strand, acceptor is 5' and donor is 3'
                                five_prime_usage = acceptor_usage
                                three_prime_usage = donor_usage

                            junction_usage[j_id] = {
                                'five_prime_usage': five_prime_usage,
                                'three_prime_usage': three_prime_usage,
                                'donor_usage': donor_usage,
                                'acceptor_usage': acceptor_usage,
                                'donor_total_reads': site_total_reads[j_data['donor']],
                                'acceptor_total_reads': site_total_reads[j_data['acceptor']]
                            }

                            # If either the donor or acceptor usage is too low, filter out the junction
                            if donor_usage < min_splice_site_usage or acceptor_usage < min_splice_site_usage:
                                filtered_out.add(j_id)
                                stats['filtered_junctions'].append({
                                    'gene_id': gene_id,
                                    'junction_id': j_id,
                                    'five_prime_usage': five_prime_usage,
                                    'three_prime_usage': three_prime_usage,
                                    'donor_usage': donor_usage,
                                    'acceptor_usage': acceptor_usage,
                                    'donor_total_reads': site_total_reads[j_data['donor']],
                                    'acceptor_total_reads': site_total_reads[j_data['acceptor']],
                                    'junction_reads': j_data['score']
                                })

                        # Remove filtered junctions
                        filtered_connected = connected - filtered_out
                        filtered_junctions += len(filtered_out)

                        # Check connectivity and split into separate ATSEs if needed
                        if len(filtered_connected) >= 2:
                            # Create a temporary subgraph with only the remaining junctions
                            temp_edges = []
                            for j_id in filtered_connected:
                                for u, v, data in G.edges(data=True):
                                    if data['junction_id'] == j_id:
                                        temp_edges.append((u, v, data))

                            sub_G = nx.Graph()
                            for u, v, data in temp_edges:
                                sub_G.add_edge(u, v, **data)

                            # Find all connected components in the subgraph
                            components = list(nx.connected_components(sub_G))

                            # Track if this ATSE was split
                            if len(components) > 1:
                                stats['split_atses'] += 1

                            # Process each connected component separately
                            for component in components:
                                # Find junctions in this component
                                component_junctions = []
                                for u, v, data in sub_G.edges(data=True):
                                    if u in component or v in component:
                                        component_junctions.append(data['junction_id'])

                                # Only create ATSE if component has ≥2 junctions
                                if len(component_junctions) >= 2:
                                    # Use temporary IDs for now
                                    temp_event_id = f"TEMP_ATSE_{event_counter}"

                                    # Calculate splice sites for this component
                                    component_splice_sites = set()
                                    for j_id in component_junctions:
                                        component_splice_sites.add(junction_data[j_id]['donor'])
                                        component_splice_sites.add(junction_data[j_id]['acceptor'])

                                    stats['junction_counts'][len(component_junctions)] += 1
                                    stats['junctions_in_atses'].update(component_junctions)
                                    stats['component_breakdown'][len(component_junctions)] += 1

                                    temp_atse_groups[temp_event_id] = {
                                        'gene_id': gene_id,
                                        'junction_ids': list(component_junctions),
                                        'num_junctions': len(component_junctions),
                                        'splice_sites': list(component_splice_sites),
                                        'filtered_junctions': list(filtered_out),
                                        'junction_usage': {j_id: junction_usage[j_id] for j_id in component_junctions},
                                        'original_group_size': len(connected),
                                        'was_split': len(components) > 1
                                    }
                                    event_counter += 1
                                else:
                                    # Component has only 1 junction - add as singleton
                                    for j_id in component_junctions:
                                        singleton_info = self.analyze_singleton_junctions(G, j_id)
                                        singleton_info['gene_id'] = gene_id
                                        stats['singleton_junctions'].append(singleton_info)
                        else:
                            # Not enough junctions after filtering - convert to singletons
                            for j_id in filtered_connected:
                                singleton_info = self.analyze_singleton_junctions(G, j_id)
                                singleton_info['gene_id'] = gene_id
                                stats['singleton_junctions'].append(singleton_info)

                        # Mark all original junctions as visited
                        visited.update(connected)
                    else:
                        # This is a singleton junction
                        singleton_info = self.analyze_singleton_junctions(G, junction_id)
                        singleton_info['gene_id'] = gene_id
                        stats['singleton_junctions'].append(singleton_info)
                        visited.add(junction_id)

        # NOW CLASSIFY, RENAME, AND REORGANIZE ATSEs
        print("Classifying and renaming ATSEs...")
        final_atse_groups, event_counts = self.classify_events(graphs, temp_atse_groups)

        # Calculate final statistics
        total_atses = len(final_atse_groups)
        junctions_used = len(stats['junctions_in_atses'])
        singleton_count = len(stats['singleton_junctions'])

        # Write singleton information to file
        with open('singleton_junctions.tsv', 'w') as f:
            # Write header
            f.write('gene_id\tjunction_id\tchromosome\tdonor_site\tacceptor_site\tstrand\tread_count\n')
            # Write data
            for j in stats['singleton_junctions']:
                f.write(f"{j['gene_id']}\t{j['junction_id']}\t{j['chromosome']}\t{j['donor_site']}\t"
                       f"{j['acceptor_site']}\t{j['strand']}\t{j['score']}\n")

        # Write filtered junction information to file
        with open('filtered_low_usage_junctions.tsv', 'w') as f:
            # Write header
            f.write('gene_id\tjunction_id\tjunction_reads\tdonor_total_reads\tdonor_usage\tacceptor_total_reads\tacceptor_usage\t5prime_usage\t3prime_usage\n')
            # Write data
            for j in stats['filtered_junctions']:
                f.write(f"{j['gene_id']}\t{j['junction_id']}\t{j['junction_reads']}\t"
                       f"{j['donor_total_reads']}\t{j['donor_usage']:.4f}\t"
                       f"{j['acceptor_total_reads']}\t{j['acceptor_usage']:.4f}\t"
                       f"{j['five_prime_usage']:.4f}\t{j['three_prime_usage']:.4f}\n")

        # Print comprehensive summary
        print(f"""
        ATSE Analysis Summary:
        ---------------------
        Total junctions in dataset: {stats['total_junctions']}
        Analyzable junctions (in genes with ≥2 junctions): {stats['analyzable_junctions']}
        Junctions filtered due to low splice site usage (<{min_splice_site_usage*100:.1f}%): {filtered_junctions}
        Junctions included in ATSEs: {junctions_used}
        Singleton junctions: {singleton_count}
        Total ATSEs found: {total_atses}

        ATSE Splitting Statistics:
        -------------------------
        Original ATSEs that were split: {stats['split_atses']}

        ATSE Size Distribution:
        ----------------------""")

        sorted_counts = dict(sorted(stats['junction_counts'].items()))
        for num_junctions, count in sorted_counts.items():
            print(f"ATSEs with {num_junctions} junctions: {count}")

        print(f"\nEvent Type Distribution:")
        print(f"  Alternative 3' splice site: {event_counts['alternative_3_prime']}")
        print(f"  Alternative 5' splice site: {event_counts['alternative_5_prime']}")  
        print(f"  Exon skipping: {event_counts['exon_skip']}")
        print(f"  Complex (multiple sites): {event_counts['complex']}")

        print(f"\nComponent size breakdown after splitting:")
        for size, count in sorted(stats['component_breakdown'].items()):
            print(f"Components with {size} junctions: {count}")

        # Generate and print gene summary
        gene_summaries = self.summarize_gene_atses(final_atse_groups)

        print(f"\nGene ATSE Summary:")
        print(f"Genes with ATSEs: {len(gene_summaries)}")

        # Show distribution of ATSEs per gene
        atses_per_gene = [summary['total_atses'] for summary in gene_summaries.values()]
        if atses_per_gene:
            print(f"Average ATSEs per gene: {sum(atses_per_gene)/len(atses_per_gene):.1f}")
            print(f"Max ATSEs in a gene: {max(atses_per_gene)}")
            print(f"Genes with only 1 ATSE: {atses_per_gene.count(1)}")
            print(f"Genes with >5 ATSEs: {len([x for x in atses_per_gene if x > 5])}")

        print(f"\nDetails of filtered junctions saved to 'filtered_low_usage_junctions.tsv'")

        return final_atse_groups, event_counts


    def classify_events(self, graphs: Dict[str, nx.Graph], atse_groups: Dict[str, Dict]):
        """
        Classify ATSE events based on splice patterns detectable from split reads.

        Event types we can reliably detect:
        - Alternative 5' splice site
        - Alternative 3' splice site  
        - Exon skipping
        - Complex (multiple donor and acceptor sites)
        """
        # Initialize counter for event types
        event_counts = {
            'alternative_3_prime': 0,
            'alternative_5_prime': 0,
            'exon_skip': 0,
            'complex': 0
        }

        # Group ATSEs by gene for ordering
        atses_by_gene = defaultdict(list)

        for event_id, group in atse_groups.items():
            gene_id = group['gene_id']
            G = graphs[gene_id]

            # Get ATSE range and strand information
            atse_start, atse_end, chromosome = self.get_atse_genomic_range(G, group['junction_ids'])

            # Get strand from first junction
            strand = None
            for j_id in group['junction_ids']:
                for _, _, data in G.edges(data=True):
                    if data['junction_id'] == j_id:
                        strand = data['strand']
                        break
                if strand:
                    break
                
            # Add genomic position info
            group['chromosome'] = chromosome
            group['atse_start'] = atse_start
            group['atse_end'] = atse_end
            group['strand'] = strand
            group['atse_length'] = atse_end - atse_start

            # Get strand-aware start position for ordering
            group['strand_aware_start'] = self.get_atse_start_position(G, group['junction_ids'], strand)

            # Count unique donor and acceptor sites
            donor_sites = len([s for s in group['splice_sites'] if s[2] == 'donor'])
            acceptor_sites = len([s for s in group['splice_sites'] if s[2] == 'acceptor'])
            num_junctions = len(group['junction_ids'])

            # Classify the event type
            if donor_sites == 1 and acceptor_sites > 1:
                group['event_type'] = 'alternative_3_prime'
            elif donor_sites > 1 and acceptor_sites == 1:
                group['event_type'] = 'alternative_5_prime'
            elif num_junctions == 3:
                # Check for exon skipping pattern
                is_exon_skip = self._check_exon_skipping(G, group['junction_ids'], strand)
                group['event_type'] = 'exon_skip' if is_exon_skip else 'complex'
            else:
                group['event_type'] = 'complex'

            # Update counter
            event_counts[group['event_type']] += 1

            # Group by gene for renaming
            atses_by_gene[gene_id].append((group['strand_aware_start'], event_id, group))

        # Reorganize ATSEs with new naming scheme
        reorganized_atses = {}

        for gene_id, gene_atses in atses_by_gene.items():
            # Sort ATSEs within gene by strand-aware position
            strand = gene_atses[0][2]['strand'] if gene_atses else '+'

            if strand == '+':
                # Positive strand: sort by increasing coordinate
                gene_atses.sort(key=lambda x: x[0])
            else:
                # Negative strand: sort by decreasing coordinate
                gene_atses.sort(key=lambda x: x[0], reverse=True)

            # Rename ATSEs within gene
            for i, (_, old_event_id, group) in enumerate(gene_atses, 1):
                new_event_id = f"{gene_id}_atse_{i}"

                # Add relative position information
                group['atse_number'] = i
                group['total_atses_in_gene'] = len(gene_atses)

                # Add distance to neighboring ATSEs
                if i > 1:
                    prev_atse = gene_atses[i-2][2]
                    group['distance_to_previous'] = abs(group['strand_aware_start'] - prev_atse['strand_aware_start'])
                else:
                    group['distance_to_previous'] = None

                if i < len(gene_atses):
                    next_atse = gene_atses[i][2]
                    group['distance_to_next'] = abs(group['strand_aware_start'] - next_atse['strand_aware_start'])
                else:
                    group['distance_to_next'] = None

                reorganized_atses[new_event_id] = group

        return reorganized_atses, event_counts

    def _check_exon_skipping(self, G: nx.Graph, junction_ids: List[str], strand: str) -> bool:
        """
        Check if three junctions form an exon skipping pattern.

        Pattern: Junction A (start-end), Junction B (start-middle), Junction C (middle-end)
        Where "middle" represents the boundaries of the skipped exon.
        """
        if len(junction_ids) != 3:
            return False

        # Get junction coordinates
        junc_coords = []
        for j_id in junction_ids:
            for u, v, data in G.edges(data=True):
                if data['junction_id'] == j_id:
                    coord1, coord2 = u[1], v[1]
                    if strand == '+':
                        start, end = min(coord1, coord2), max(coord1, coord2)
                    else:
                        start, end = max(coord1, coord2), min(coord1, coord2)
                    junc_coords.append((start, end, j_id))
                    break
                
        # Sort junctions by their start coordinate (considering strand)
        junc_coords.sort(key=lambda x: x[0] if strand == '+' else -x[0])

        # Check for exon skipping pattern
        # We need: J1(A-C), J2(A-B), J3(B-C) where B is the skipped exon
        for i in range(len(junc_coords)):
            for j in range(len(junc_coords)):
                for k in range(len(junc_coords)):
                    if i == j or j == k or i == k:
                        continue
                    
                    j1_start, j1_end, _ = junc_coords[i]
                    j2_start, j2_end, _ = junc_coords[j]
                    j3_start, j3_end, _ = junc_coords[k]

                    # Check if we have the pattern: J1 spans, J2 and J3 define boundaries
                    if strand == '+':
                        if (j1_start == j2_start and j1_end == j3_end and j2_end == j3_start):
                            return True
                    else:
                        if (j1_end == j2_end and j1_start == j3_start and j2_start == j3_end):
                            return True

        return False

    def get_atse_genomic_range(self, graph: nx.Graph, junction_ids: List[str]) -> Tuple[int, int, str]:
        """
        Get the genomic range (start, end) of an ATSE and its chromosome.

        Returns:
            Tuple of (start, end, chromosome) where start < end regardless of strand
        """
        coordinates = []
        chromosome = None

        for j_id in junction_ids:
            for u, v, data in graph.edges(data=True):
                if data['junction_id'] == j_id:
                    coord1, coord2 = u[1], v[1]
                    chromosome = u[0]
                    coordinates.extend([coord1, coord2])
                    break
                
        if not coordinates:
            return 0, 0, ""

        return min(coordinates), max(coordinates), chromosome

    def get_atse_start_position(self, graph: nx.Graph, junction_ids: List[str], strand: str) -> int:
        """
        Get the strand-aware start position of an ATSE.

        For positive strand: returns leftmost coordinate
        For negative strand: returns rightmost coordinate
        """
        coordinates = []

        for j_id in junction_ids:
            for u, v, data in graph.edges(data=True):
                if data['junction_id'] == j_id:
                    coordinates.extend([u[1], v[1]])
                    break
                
        if not coordinates:
            return 0

        if strand == '+':
            return min(coordinates)
        else:
            return max(coordinates)


    def summarize_gene_atses(self, atse_groups: Dict[str, Dict]) -> Dict[str, Dict]:
        """
        Create a summary of ATSEs per gene with their relative positions.
        """
        gene_summaries = defaultdict(lambda: {
            'total_atses': 0,
            'event_types': defaultdict(int),
            'atses': [],
            'gene_atse_range': None
        })

        for event_id, group in atse_groups.items():
            gene_id = group['gene_id']
            summary = gene_summaries[gene_id]

            summary['total_atses'] += 1
            summary['event_types'][group['event_type']] += 1
            summary['atses'].append({
                'event_id': event_id,
                'event_number': group.get('atse_number', 0),
                'start': group.get('atse_start', 0),
                'end': group.get('atse_end', 0),
                'event_type': group['event_type'],
                'num_junctions': group['num_junctions']
            })

            # Update gene range
            if summary['gene_atse_range'] is None:
                summary['gene_atse_range'] = [group.get('atse_start', 0), group.get('atse_end', 0)]
            else:
                summary['gene_atse_range'][0] = min(summary['gene_atse_range'][0], group.get('atse_start', 0))
                summary['gene_atse_range'][1] = max(summary['gene_atse_range'][1], group.get('atse_end', 0))

        # Convert defaultdict to regular dict and calculate span
        for gene_id, summary in gene_summaries.items():
            if summary['gene_atse_range']:
                summary['gene_atse_span'] = summary['gene_atse_range'][1] - summary['gene_atse_range'][0]
        return dict(gene_summaries)

    def save_atse_file(self, atse_groups: Dict[str, Dict], junctions: Dict[str, Dict], file_name: str):
        """
        Save ATSE groups to a tab-delimited file with gzip compression, including positional information.
        """
        import gzip
        from collections import defaultdict

        required_fields = {'gene_id', 'num_junctions', 'event_type', 'junction_ids'}

        # Ensure file has .gz extension
        if not file_name.endswith('.gz'):
            file_name = file_name + '.gz'

        # Define column groups with positional information
        event_columns = [
            "event_id", "gene_id", "gene_name", "gene_types",
            "num_junctions", "event_type", "chromosome", "strand",
            "atse_start", "atse_end", "atse_length",
            "atse_number", "total_atses_in_gene",
            "distance_to_previous", "distance_to_next"
        ]

        transcript_columns = [
            "transcripts", "both_ends_transcripts", "only_5_prime_transcripts", 
            "only_3_prime_transcripts", "transcript_types", "annotation_status",
            "perfect_match_5_prime", "perfect_match_3_prime"
        ]

        junction_columns = [
            "junction_id", "chrom", "start", "end", "strand", "cells", "total_score"
        ]

        usage_columns = [
            "five_prime_usage", "three_prime_usage", "donor_usage", "acceptor_usage",
            "donor_total_reads", "acceptor_total_reads"
        ]

        sequence_columns = [
            "splice_motif", "donor_seq", "acceptor_seq"
        ]

        position_columns = [
            "position_off_5_prime", "position_off_3_prime"
        ]

        # Combine all columns
        all_columns = event_columns + transcript_columns + junction_columns + usage_columns + sequence_columns + position_columns

        # Track statistics
        stats = defaultdict(int)

        try:
            with gzip.open(file_name, 'wt') as f:
                # Write header
                f.write("\t".join(all_columns) + "\n")

                # Write data
                for event_id, group in atse_groups.items():
                    # Verify all required fields are present
                    missing_fields = required_fields - set(group.keys())
                    if missing_fields:
                        print(f"Warning: Event {event_id} missing required fields: {missing_fields}")
                        stats["missing_fields"] += 1
                        continue
                    
                    # Get junction usage data if available
                    junction_usage = group.get('junction_usage', {})

                    try:
                        # For each junction in the ATSE
                        for junction_id in group['junction_ids']:
                            if junction_id not in junctions:
                                print(f"Warning: Junction {junction_id} not found in annotations")
                                stats["missing_junctions"] += 1
                                continue
                            
                            j_data = junctions[junction_id]
                            stats["junctions_written"] += 1

                            # Prepare row data
                            row_data = {}

                            # Event data with positional information
                            row_data["event_id"] = event_id
                            row_data["gene_id"] = group['gene_id']

                            # Handle gene names
                            gene_names = j_data.get('gene_names', [])
                            row_data["gene_name"] = '|'.join(str(name) for name in gene_names) if gene_names else 'NA'

                            # Handle gene types
                            gene_types = j_data.get('gene_types', [])
                            row_data["gene_types"] = '|'.join(str(gtype) for gtype in gene_types) if gene_types else 'NA'

                            row_data["num_junctions"] = group['num_junctions']
                            row_data["event_type"] = group['event_type']

                            # Add positional information
                            row_data["chromosome"] = group.get('chromosome', 'NA')
                            row_data["strand"] = group.get('strand', 'NA')
                            row_data["atse_start"] = group.get('atse_start', 'NA')
                            row_data["atse_end"] = group.get('atse_end', 'NA')
                            row_data["atse_length"] = group.get('atse_length', 'NA')
                            row_data["atse_number"] = group.get('atse_number', 'NA')
                            row_data["total_atses_in_gene"] = group.get('total_atses_in_gene', 'NA')
                            row_data["distance_to_previous"] = group.get('distance_to_previous', 'NA')
                            row_data["distance_to_next"] = group.get('distance_to_next', 'NA')

                            # Transcript data
                            for field in ["transcripts", "both_ends_transcripts", "only_5_prime_transcripts", 
                                         "only_3_prime_transcripts", "perfect_match_5_prime", "perfect_match_3_prime"]:
                                values = j_data.get(field, [])
                                row_data[field] = ','.join(str(t) for t in values) if values else 'NA'

                            # Transcript types
                            transcript_types = j_data.get('transcript_types', [])
                            row_data["transcript_types"] = ','.join(str(t) for t in transcript_types) if transcript_types else 'NA'

                            row_data["annotation_status"] = j_data.get('annotation_status', 'NA')

                            # Junction data
                            row_data["junction_id"] = junction_id
                            for field in ["chrom", "start", "end", "strand", "cells", "total_score"]:
                                row_data[field] = j_data.get(field, 'NA')

                            # Usage data
                            usage_data = junction_usage.get(junction_id, {})
                            for field in ["five_prime_usage", "three_prime_usage", "donor_usage", "acceptor_usage"]:
                                value = usage_data.get(field, 'NA')
                                row_data[field] = f"{value:.4f}" if isinstance(value, float) else 'NA'

                            # Total reads data
                            for field in ["donor_total_reads", "acceptor_total_reads"]:
                                value = usage_data.get(field, 'NA')
                                row_data[field] = f"{value}" if isinstance(value, (int, float)) else 'NA'

                            # Sequence data
                            for field in ["splice_motif", "donor_seq", "acceptor_seq"]:
                                row_data[field] = j_data.get(field, 'NA')

                            # Position data
                            for field in ["position_off_5_prime", "position_off_3_prime"]:
                                row_data[field] = j_data.get(field, 'NA')

                            # Write the row
                            f.write("\t".join(str(row_data.get(col, 'NA')) for col in all_columns) + "\n")

                    except Exception as e:
                        print(f"Warning: Error writing event {event_id}: {str(e)}")
                        stats["write_errors"] += 1
                        continue
                    
            # Print statistics
            print(f"\nATSEs successfully saved to {file_name}")
            print(f"Wrote {len(atse_groups)} ATSE groups with {stats['junctions_written']} junctions")
            if stats["missing_fields"] > 0:
                print(f"Skipped {stats['missing_fields']} events due to missing fields")
            if stats["missing_junctions"] > 0:
                print(f"Skipped {stats['missing_junctions']} junctions not found in annotations")
            if stats["write_errors"] > 0:
                print(f"Encountered {stats['write_errors']} errors while writing events")

        except IOError as e:
            print(f"Error: Could not write to file {file_name}: {str(e)}")
            raise
        except Exception as e:
            print(f"Error: Unexpected error while saving ATSEs: {str(e)}")
            raise