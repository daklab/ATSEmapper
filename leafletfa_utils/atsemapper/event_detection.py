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
            gene_id = j_data['gene_ids'][0]
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
        
    def sample_unincluded_junctions(self, graphs: Dict[str, nx.Graph], included_junctions: set, sample_size: int = 500) -> List[Dict]:
        """
        Sample unincluded junctions and verify they don't share splice sites.
        Returns list of any problematic cases found.
        """
        all_unincluded = []
        for gene_id, G in graphs.items():
            gene_junctions = {data['junction_id'] for _, _, data in G.edges(data=True)}
            unincluded = gene_junctions - included_junctions
            all_unincluded.extend((gene_id, j_id) for j_id in unincluded)

        # Sample junctions
        if len(all_unincluded) > sample_size:
            sampled = random.sample(all_unincluded, sample_size)
        else:
            sampled = all_unincluded

        problematic = []

        # Check each sampled junction
        for gene_id, junction_id in sampled:
            G = graphs[gene_id]

            # Get junction's splice sites
            current_donor = None
            current_acceptor = None
            for u, v, data in G.edges(data=True):
                if data['junction_id'] == junction_id:
                    current_donor = u[1]  # Just the coordinate
                    current_acceptor = v[1]
                    break
                
            if current_donor is None or current_acceptor is None:
                print(f"Warning: Could not find coordinates for {junction_id}")
                continue
            
            # Check all other junctions in the same gene
            for u, v, data in G.edges(data=True):
                other_id = data['junction_id']
                if other_id != junction_id:
                    other_donor = u[1]
                    other_acceptor = v[1]

                    # Check for exact coordinate matches
                    if (current_donor == other_donor or 
                        current_donor == other_acceptor or
                        current_acceptor == other_donor or 
                        current_acceptor == other_acceptor):

                        # For debugging, print the exact match found
                        problematic.append({
                            'gene_id': gene_id,
                            'junction_id': junction_id,
                            'junction_coords': (current_donor, current_acceptor),
                            'shares_site_with': other_id,
                            'other_coords': (other_donor, other_acceptor),
                            'shared_coord': current_donor if (current_donor == other_donor or current_donor == other_acceptor) else current_acceptor
                        })

        if problematic:
            print("\nDetailed analysis of problematic cases:")
            for case in problematic[:5]:
                print(f"\nJunction: {case['junction_id']} ({case['junction_coords']})")
                print(f"Shares with: {case['shares_site_with']} ({case['other_coords']})")
                print(f"Shared coordinate: {case['shared_coord']}")

        return problematic
    
    def find_atse_groups(self, graphs: Dict[str, nx.Graph], min_splice_site_usage: float = 0.01) -> Dict[str, Dict]:
        """
        Find alternative transcript splicing events (ATSEs) with splice site usage filtering.
    
        Args:
            graphs: Dictionary of gene graphs
            min_splice_site_usage: Minimum proportion of reads a junction must have at a splice site
                               compared to total reads at that site (default: 0.01 or 1%)
    
        Returns:
            Dictionary of ATSE groups and sorted counts
        """
        atse_groups = {}
        event_counter = 0
        filtered_junctions = 0  # Track filtered junctions
    
        # Statistics tracking
        stats = {
            'total_junctions': sum(G.number_of_edges() for G in graphs.values()),
            'analyzable_junctions': 0,
            'junctions_in_atses': set(),
            'junction_counts': defaultdict(int),
            'singleton_junctions': [],  # Store info about singleton junctions
            'filtered_junctions': []    # Track junctions filtered due to low splice site usage
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

                        # Only create an ATSE if at least 2 junctions remain after filtering
                        if len(filtered_connected) >= 2:
                            event_id = f"ATSE_{event_counter}"
                        
                            # Recalculate splice sites based on filtered junctions
                            filtered_splice_sites = set()
                            for j_id in filtered_connected:
                                filtered_splice_sites.add(junction_data[j_id]['donor'])
                                filtered_splice_sites.add(junction_data[j_id]['acceptor'])

                            stats['junction_counts'][len(filtered_connected)] += 1
                            stats['junctions_in_atses'].update(filtered_connected)

                            atse_groups[event_id] = {
                                'gene_id': gene_id,
                                'junction_ids': list(filtered_connected),
                                'num_junctions': len(filtered_connected),
                                'splice_sites': list(filtered_splice_sites),
                                'filtered_junctions': list(filtered_out),
                                'junction_usage': {j_id: junction_usage[j_id] for j_id in filtered_connected}
                            }
                            event_counter += 1
                        else:
                            # Add the remaining junctions as singletons if they don't form an ATSE anymore
                            for j_id in filtered_connected:
                                singleton_info = self.analyze_singleton_junctions(G, j_id)
                                singleton_info['gene_id'] = gene_id
                                stats['singleton_junctions'].append(singleton_info)

                        # Mark all junctions as visited
                        visited.update(connected)
                    else:
                        # This is a singleton junction
                        singleton_info = self.analyze_singleton_junctions(G, junction_id)
                        singleton_info['gene_id'] = gene_id
                        stats['singleton_junctions'].append(singleton_info)
                        visited.add(junction_id)

        # Calculate final statistics
        total_atses = len(atse_groups)
        junctions_used = len(stats['junctions_in_atses'])
        singleton_count = len(stats['singleton_junctions'])
    
        # Perform sanity check on sample of unincluded junctions
        problematic = self.sample_unincluded_junctions(graphs, stats['junctions_in_atses'])

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

        print(f"""
                ATSE Analysis Summary:
                ---------------------
                Total junctions in dataset: {stats['total_junctions']}
                Analyzable junctions (in genes with ≥2 junctions): {stats['analyzable_junctions']}
                Junctions filtered due to low splice site usage (<{min_splice_site_usage*100:.1f}%): {filtered_junctions}
                Junctions included in ATSEs: {junctions_used}
                Singleton junctions: {singleton_count}
                Total ATSEs found: {total_atses}

                ATSE Size Distribution:
                ----------------------""")

        sorted_counts = dict(sorted(stats['junction_counts'].items()))
        for num_junctions, count in sorted_counts.items():
            print(f"ATSEs with {num_junctions} junctions: {count}")

        if problematic:
            print(f"\nWARNING: Found {len(problematic)} potentially problematic cases in random sampling")
            print("First few examples:")
            for case in problematic[:5]:
                print(f"Junction {case['junction_id']} in gene {case['gene_id']} shares site with {case['shares_site_with']}")
        else:
            print("\nRandom sampling verification: OK - no shared splice sites found in sampled junctions")

        print(f"\nDetails of filtered junctions saved to 'filtered_low_usage_junctions.tsv'")

        return atse_groups, sorted_counts
        
    def classify_events(self, graphs: Dict[str, nx.Graph], atse_groups: Dict[str, Dict]):
        # Initialize counter for event types
        event_counts = {
            'alternative_3_prime': 0,
            'alternative_5_prime': 0,
            'exon_skip': 0,
            'complex': 0
        }

        for event_id, group in atse_groups.items():
            G = graphs[group['gene_id']]

            # Count unique donor and acceptor sites
            donor_sites = len([s for s in group['splice_sites'] if s[2] == 'donor'])
            acceptor_sites = len([s for s in group['splice_sites'] if s[2] == 'acceptor'])

            # Need a window on how far the alternative splice sites are 
            if donor_sites == 1 and acceptor_sites > 1:
                group['event_type'] = 'alternative_3_prime'
            elif donor_sites > 1 and acceptor_sites == 1:
                group['event_type'] = 'alternative_5_prime'
            elif len(group['junction_ids']) == 3:
                # Get junction coordinates and strand
                junc_coords = []
                strand = None
                for j_id in group['junction_ids']:
                    for u, v, data in G.edges(data=True):
                        if data['junction_id'] == j_id:
                            coord1, coord2 = u[1], v[1]
                            strand = u[2]  # Get strand from node tuple (assuming format: (gene_id, position, strand))
                            junc_coords.append((coord1, coord2))
                            break
                        
                # Sort junctions based on start position, accounting for strand
                if strand == '+':
                    # For positive strand, smaller coordinate is start
                    junc_coords.sort(key=lambda x: min(x))
                else:
                    # For negative strand, larger coordinate is start
                    junc_coords.sort(key=lambda x: -max(x))

                # Check for exon skipping by comparing starts and ends based on strand
                is_exon_skip = False
                if strand == '+':
                    # Positive strand: compare smallest coordinates for starts, largest for ends
                    if (min(junc_coords[0]) == min(junc_coords[1]) and  # J1 start == J2 start
                        max(junc_coords[1]) == max(junc_coords[2])):    # J2 end == J3 end
                        is_exon_skip = True
                else:
                    # Negative strand: compare largest coordinates for starts, smallest for ends
                    if (max(junc_coords[0]) == max(junc_coords[1]) and  # J1 start == J2 start
                        min(junc_coords[1]) == min(junc_coords[2])):    # J2 end == J3 end
                        is_exon_skip = True

                group['event_type'] = 'exon_skip' if is_exon_skip else 'complex'
            else:
                group['event_type'] = 'complex'

            # Update counter
            event_counts[group['event_type']] += 1

        return atse_groups, event_counts
    
    def save_atse_file(self, atse_groups: Dict[str, Dict], junctions: Dict[str, Dict], file_name: str):
        """
        Save ATSE groups to a tab-delimited file with gzip compression, including junction annotations.

        Args:
            atse_groups: Dictionary of ATSE events
            junctions: Dictionary of junction annotations
            file_name: Output file path (will append .gz if not present)
        """
        import gzip

        required_fields = {'gene_id', 'num_junctions', 'event_type', 'junction_ids'}

        # Ensure file has .gz extension
        if not file_name.endswith('.gz'):
            file_name = file_name + '.gz'

        try:
            with gzip.open(file_name, 'wt') as f:  # 'wt' for write text mode
                # Write header
                f.write("event_id\tgene_id\tgene_name\tgene_types\t"
                   "transcripts\tboth_ends_transcripts\tonly_5_prime_transcripts\tonly_3_prime_transcripts\t"
                   "transcript_types\tnum_junctions\tevent_type\tannotation_status\t"
                   "junction_id\tchrom\tstart\tend\tstrand\tcells\ttotal_score\t"
                   "five_prime_usage\tthree_prime_usage\tdonor_usage\tacceptor_usage\t"
                   "donor_total_reads\tacceptor_total_reads\t"  # New columns
                   "splice_motif\tdonor_seq\tacceptor_seq\t"
                   "position_off_5_prime\tposition_off_3_prime\n")

                # Write data
                for event_id, group in atse_groups.items():
                    # Verify all required fields are present
                    missing_fields = required_fields - set(group.keys())
                    if missing_fields:
                        print(f"Warning: Event {event_id} missing required fields: {missing_fields}")
                        continue
                    
                    # Get junction usage data if available
                    junction_usage = group.get('junction_usage', {})

                    try:
                        # For each junction in the ATSE
                        for junction_id in group['junction_ids']:
                            if junction_id not in junctions:
                                print(f"Warning: Junction {junction_id} not found in annotations")
                                continue

                            j_data = junctions[junction_id]

                            # Get usage values for this junction
                            usage_data = junction_usage.get(junction_id, {})
                            five_prime_usage = usage_data.get('five_prime_usage', 'NA')
                            three_prime_usage = usage_data.get('three_prime_usage', 'NA')
                            donor_usage = usage_data.get('donor_usage', 'NA')
                            acceptor_usage = usage_data.get('acceptor_usage', 'NA')

                            # Get donor and acceptor total reads
                            donor_total_reads = usage_data.get('donor_total_reads', 'NA')
                            acceptor_total_reads = usage_data.get('acceptor_total_reads', 'NA')

                            # Format usage values
                            five_prime_usage_str = f"{five_prime_usage:.4f}" if isinstance(five_prime_usage, float) else 'NA'
                            three_prime_usage_str = f"{three_prime_usage:.4f}" if isinstance(three_prime_usage, float) else 'NA'
                            donor_usage_str = f"{donor_usage:.4f}" if isinstance(donor_usage, float) else 'NA'
                            acceptor_usage_str = f"{acceptor_usage:.4f}" if isinstance(acceptor_usage, float) else 'NA'

                            # Format total reads values
                            donor_total_reads_str = f"{donor_total_reads}" if isinstance(donor_total_reads, (int, float)) else 'NA'
                            acceptor_total_reads_str = f"{acceptor_total_reads}" if isinstance(acceptor_total_reads, (int, float)) else 'NA'

                            # Handle gene names - join with pipe if multiple names exist
                            gene_names = j_data.get('gene_names', [])
                            gene_names_str = '|'.join(str(name) for name in gene_names) if gene_names else 'NA'

                            # Handle gene types - new column
                            gene_types = j_data.get('gene_types', [])
                            gene_types_str = '|'.join(str(gtype) for gtype in gene_types) if gene_types else 'NA'

                            # Handle transcripts - join with comma or return NA if empty
                            transcripts = j_data.get('transcripts', [])
                            transcripts_str = ','.join(str(t) for t in transcripts) if transcripts else 'NA'

                            # Handle the three transcript categories - new columns
                            both_ends = j_data.get('both_ends_transcripts', [])
                            both_ends_str = ','.join(str(t) for t in both_ends) if both_ends else 'NA'

                            only_5_prime = j_data.get('only_5_prime_transcripts', [])
                            only_5_prime_str = ','.join(str(t) for t in only_5_prime) if only_5_prime else 'NA'

                            only_3_prime = j_data.get('only_3_prime_transcripts', [])
                            only_3_prime_str = ','.join(str(t) for t in only_3_prime) if only_3_prime else 'NA'

                            # Handle transcripts types - join with comma or return NA if empty
                            transcript_types = j_data.get('transcript_types', [])
                            transcript_types_str = ','.join(str(t) for t in transcript_types) if transcript_types else 'NA'

                            # Write line with all information, including new columns
                            f.write(f"{event_id}\t"
                               f"{group['gene_id']}\t"
                               f"{gene_names_str}\t"
                               f"{gene_types_str}\t"
                               f"{transcripts_str}\t"
                               f"{both_ends_str}\t"
                               f"{only_5_prime_str}\t"
                               f"{only_3_prime_str}\t"
                               f"{transcript_types_str}\t"
                               f"{group['num_junctions']}\t"
                               f"{group['event_type']}\t"
                               f"{j_data.get('annotation_status', 'NA')}\t"
                               f"{junction_id}\t"
                               f"{j_data.get('chrom', 'NA')}\t"
                               f"{j_data.get('start', 'NA')}\t"
                               f"{j_data.get('end', 'NA')}\t"
                               f"{j_data.get('strand', 'NA')}\t"
                               f"{j_data.get('cells', 'NA')}\t"
                               f"{j_data.get('total_score', 'NA')}\t"
                               f"{five_prime_usage_str}\t"
                               f"{three_prime_usage_str}\t"
                               f"{donor_usage_str}\t"
                               f"{acceptor_usage_str}\t"
                               f"{donor_total_reads_str}\t"  # New column
                               f"{acceptor_total_reads_str}\t"  # New column
                               f"{j_data.get('splice_motif', 'NA')}\t"
                               f"{j_data.get('donor_seq', 'NA')}\t"
                               f"{j_data.get('acceptor_seq', 'NA')}\t"
                               f"{j_data.get('position_off_5_prime', 'NA')}\t"
                               f"{j_data.get('position_off_3_prime', 'NA')}\n"
                            )
                    except Exception as e:
                        print(f"Warning: Error writing event {event_id}: {str(e)}")
                        continue
                    
            print(f"ATSEs successfully saved to {file_name}")
            print(f"Wrote {len(atse_groups)} ATSE events")
        
        except IOError as e:
            print(f"Error: Could not write to file {file_name}: {str(e)}")
            raise
        except Exception as e:
            print(f"Error: Unexpected error while saving ATSEs: {str(e)}")
            raise