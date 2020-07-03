from PW_explorer.load_worlds import load_worlds
from PW_explorer.run_clingo import run_clingo
from PW_explorer.visualize import PWEVisualization
from PW_explorer.helper import pw_slicer, rel_slicer, rel_name_remapper
from PW_explorer.export import PWEExport
import pandas as pd

import networkx as nx
from nxpd import draw
import nxpd

from IPython.display import HTML as html_print

def cstr(s, color='black'):
    return "<text style=color:{}>{}</text>".format(color, s)

def rule_to_node(atom, occ):
    return "{}[{}]".format(atom.strip('"'), occ)

def head_to_node(head):
    return remove_double_quotes(head)

def remove_double_quotes(r):
    return r.strip('"')

def print_colored_query_string(pw_rel_dfs):
    
    heads = {h: {} for h in set(pw_rel_dfs['ruleH_1']['HEAD'])}
    rules = {}
    for _, p in pw_rel_dfs['ruleOcc_2'][['ATOM', 'OCC']].iterrows():
        rules[(p['ATOM'], int(p['OCC']))] = {}
    
    true_heads = set(pw_rel_dfs['ruleHTrue_1']['HEAD'])
    true_rules = set([])
    for _, p in pw_rel_dfs['ruleOccTrue_2'][['ATOM', 'OCC']].iterrows():
        true_rules.add((p['ATOM'], int(p['OCC'])))
    
    for _, row in pw_rel_dfs['newHArc_3'].iterrows():
        pos = int(row['POS'])
        h = row['HEAD']
        var = row['VAR']
        heads[h][pos] = var
    for _, row in pw_rel_dfs['newArc_4'].iterrows():
        rule = (row['ATOM'], int(row['OCC']))
        pos = int(row['POS'])
        var = row['VAR']
        rules[rule][pos] = var
        
#     print(heads, rules, true_heads, true_rules)
    
    heads_htmls = []
    for h in sorted(heads.keys()):
        variables = [remove_double_quotes(heads[h][k]) for k in sorted(heads[h].keys())]
        heads_htmls.append(cstr('{}({})'.format(head_to_node(h), ','.join(variables)), 'green' if h in true_heads else 'red'))
    
    
    rules_htmls = []
    for (r, occ) in sorted(rules.keys()):
        variables = [remove_double_quotes(rules[(r,occ)][k]) for k in sorted(rules[(r,occ)].keys())]
        rules_htmls.append(cstr('{}({})'.format(head_to_node(r), ','.join(variables)), 'green' if (r,occ) in true_rules else 'red'))
    
#     display(html_print(cstr(' '.join(heads_htmls + rules_htmls))))
    
    
    heads_htmls_combined = cstr(' ; '.join(heads_htmls))
    rules_htmls_combined = cstr(', '.join(rules_htmls))
    display(html_print(cstr('{} :- {}.'.format(heads_htmls_combined, rules_htmls_combined))))
    
    eqls = []
    for _, row in pw_rel_dfs['eqOrd_3'].iterrows():
        eqls.append('{} = {}'.format(row['VAR1'], row['VAR2']))
    print('where' if len(eqls) > 0 else '')
    print(' and '.join(eqls))
    
    
def get_pattern_graph(pw_rel_dfs):
    g = nx.MultiDiGraph()
    true_heads = set(pw_rel_dfs['ruleHTrue_1']['HEAD'])
    if len(true_heads) > 0:
        print("Success Pattern:")
    else:
        print("Failure Pattern:")
    
    if 'e_3' in pw_rel_dfs:
        for _, row in pw_rel_dfs['e_3'].iterrows():
            g.add_edge(row['NODE1'], row['NODE2'], color='darkgreen', label='{}{}'.format('e', row['OCC']))
    if 'ne_3' in pw_rel_dfs:
        for _, row in pw_rel_dfs['ne_3'].iterrows():
            g.add_edge(row['NODE1'], row['NODE2'], color='red', style='dotted', label='{}{}'.format('e', row['OCC']))
    
    return g

def get_pattern_graph2(pw_rel_dfs, no_node_labels=False, chain_eq_node_labels=True, silent=False):
    g = nx.DiGraph()
    true_heads = set(pw_rel_dfs['ruleHTrue_1']['HEAD'])
    g.graph['rankdir'] = 'LR'
    
    if not silent:
        if len(true_heads) > 0:
            print("Success Pattern:")
        else:
            print("Failure Pattern:")
    
    if 'e_3' in pw_rel_dfs:
        for _, row in pw_rel_dfs['e_3'].iterrows():
            if (row['NODE1'], row['NODE2']) in g.edges:
                g.edges[(row['NODE1'], row['NODE2'])]['label'] += ';{}{}'.format('e', row['OCC'])
            else:
                g.add_edge(row['NODE1'], row['NODE2'], color='darkgreen', label='{}{}'.format('e', row['OCC']),
                           fontname = "helvetica")
    if 'ne_3' in pw_rel_dfs:
        for _, row in pw_rel_dfs['ne_3'].iterrows():
            if (row['NODE1'], row['NODE2']) in g.edges:
                g.edges[(row['NODE1'], row['NODE2'])]['label'] += ';{}{}'.format('e', row['OCC'])
            else:
                g.add_edge(row['NODE1'], row['NODE2'], color='red', style='dotted', label='{}{}'.format('e', row['OCC']),
                           fontname = "helvetica")
    
    
    for n in g.nodes:
        if no_node_labels:
            g.nodes[n]['label'] = '  '
        g.nodes[n]['fontname'] = "helvetica"
        g.nodes[n]['style'] = "rounded"
        g.nodes[n]['shape'] = "box"
        if chain_eq_node_labels:
            eq_grps = eq_groups(pw_rel_dfs, single_grps=False)
            eq_grps = [set([s for s in eq_grp]) for eq_grp in eq_grps]
            for eq_grp in eq_grps:
                if n in eq_grp:
                    g.nodes[n]['label'] = '='.join(map(remove_double_quotes, sorted(list(eq_grp))))
    
    return g

def get_pattern_graph3(pw_rel_dfs, pw_obj, no_node_labels=False, chain_eq_node_labels=True, silent=False, q_head='none', q_head_color='blue'):
    """
    q_head options: none, hyperedge, binary-edge, binary-node-highlight
    """
    g = nx.DiGraph()
    true_heads = set(pw_rel_dfs['ruleHTrue_1']['HEAD'])
    g.graph['rankdir'] = 'LR'
    
    if not silent:
        if len(true_heads) > 0:
            print("Success Pattern:")
        else:
            print("Failure Pattern:")
    
    if 'e_3' in pw_rel_dfs:
        for _, row in pw_rel_dfs['e_3'].iterrows():
            if (row['NODE1'], row['NODE2']) in g.edges:
                g.edges[(row['NODE1'], row['NODE2'])]['label'] += ';{}{}'.format('e', row['OCC'])
            else:
                g.add_edge(row['NODE1'], row['NODE2'], color='darkgreen', label='{}{}'.format('e', row['OCC']),
                           fontname = "helvetica")
    if 'ne_3' in pw_rel_dfs:
        for _, row in pw_rel_dfs['ne_3'].iterrows():
            if (row['NODE1'], row['NODE2']) in g.edges:
                g.edges[(row['NODE1'], row['NODE2'])]['label'] += ';{}{}'.format('e', row['OCC'])
            else:
                g.add_edge(row['NODE1'], row['NODE2'], color='red', style='dotted', label='{}{}'.format('e', row['OCC']),
                           fontname = "helvetica")
    
    for n in g.nodes:
        if no_node_labels:
            g.nodes[n]['label'] = '  '
        g.nodes[n]['fontname'] = "Helvetica-Narrow"
        g.nodes[n]['style'] = "rounded"
        g.nodes[n]['shape'] = "box"
        if chain_eq_node_labels:
            eq_grps = eq_groups(pw_rel_dfs, single_grps=False)
            eq_grps = [set([s for s in eq_grp]) for eq_grp in eq_grps]
            for eq_grp in eq_grps:
                if n in eq_grp:
                    g.nodes[n]['label'] = '='.join(map(remove_double_quotes, sorted(list(eq_grp))))
    
    if q_head == 'hyperedge':
        heads, true_heads = get_query_heads_dict(pw_obj)
        for head, head_desc in heads.items():
            g.add_node(head, color='skyblue', shape='oval', fontname='helvetica', style='"filled,rounded"', fillcolor='skyblue')
            for pos, var in head_desc.items():
                if (head, var) in g.edges:
                    g.edges[(head,var)]['label'] += f';{str(pos)}'
                else:
                    g.add_edge(head, var, label=str(pos), color=q_head_color, style='dotted', constraint='false')
    elif q_head == 'none':
        pass
    elif q_head == 'binary-edge':
        heads, true_heads = get_query_heads_dict(pw_obj)
        for head, head_desc in heads.items():
            assert len(head_desc.items()) == 2
            [var1, var2] = [head_desc[pos] for pos in sorted(head_desc.keys())]
            g.add_edge(var1, var2, label=str(head), color=q_head_color, style='dotted', constraint='false')
    elif q_head == 'binary-node-highlight':
        heads, true_heads = get_query_heads_dict(pw_obj)
        for head, head_desc in heads.items():
            assert len(head_desc.items()) == 2
            [var1, var2] = [head_desc[pos] for pos in sorted(head_desc.keys())]
            if var1 == var2:
                g.nodes[var1]['shape'] = 'doubleoctagon'
                g.nodes[var1]['style'] = "rounded"
                g.nodes[var1]['color'] = q_head_color
            else:
                g.nodes[var1]['shape'] = 'doublecircle'
                g.nodes[var1]['style'] = "solid"
                g.nodes[var1]['color'] = 'blue'
                g.nodes[var2]['shape'] = 'doubleoctagon'
                g.nodes[var2]['style'] = "solid"
                g.nodes[var2]['color'] = q_head_color
    
    return g

def get_incidence_graph(pw_rel_dfs):
    
    g = nx.MultiDiGraph()
    g.graph['rankdir'] = 'LR'
    
    for _, row in pw_rel_dfs['ruleH_1'].iterrows():
        g.add_node(head_to_node(row['HEAD']), rank='min', color='red', shape='box', style='filled')
    for _, row in pw_rel_dfs['ruleHTrue_1'].iterrows():
        g.add_node(head_to_node(row['HEAD']), rank='min', color='green', style='filled')
    for _, row in pw_rel_dfs['ruleOcc_2'].iterrows():
        g.add_node(rule_to_node(row['ATOM'], row['OCC']), color='red', style='filled')
    for _, row in pw_rel_dfs['ruleOccTrue_2'].iterrows():
        g.add_node(rule_to_node(row['ATOM'], row['OCC']), color='green', style='filled')
    
    
    for _, row in pw_rel_dfs['newHArc_3'].iterrows():
        g.add_edge(row['VAR'], head_to_node(row['HEAD']), label=row['POS'])
    for _, row in pw_rel_dfs['newArc_4'].iterrows():
        g.add_edge(row['VAR'], rule_to_node(row['ATOM'], row['OCC']), label=row['POS'])
    
    for _, row in pw_rel_dfs['eqOrd_3'].iterrows():
        g.add_edge(row['VAR2'], row['VAR1'], constraint='false', color='green', style='dotted')
    
    return g

def eq_groups(pw_rel_dfs, single_grps=False):
    
    p = pw_rel_dfs['eqOrdMinimal_3'].groupby(['VAR1'])
    new_vars = set(pw_rel_dfs['newVar_2']['VAR'])
    groups = []
    for var1, group in p.groups.items():
        new_vars.remove(var1)
        groups.append([var1]+[pw_rel_dfs['eqOrdMinimal_3'].loc[i]['VAR2'] for i in group])
    if single_grps:
        for v in new_vars:
            groups.append([v])
    return groups
    

def chain_eqs_rules(pw_rel_dfs):
    
    eqls = []
    groups = eq_groups(pw_rel_dfs)
    for group in groups:
        eqls.append('='.join(map(remove_double_quotes, group)))
    
    return eqls

def eq_rules(pw_rel_dfs):
    
    eqls = []
    groups = eq_groups(pw_rel_dfs)
    
    for group in groups:
        var1 = group[0]
        for var2 in group[1:]:
            eqls.append('{}={}'.format(remove_double_quotes(var1), remove_double_quotes(var2)))
    
    return eqls

def neq_rules(pw_rel_dfs):
    neqls = []
    for _, row in pw_rel_dfs['neqOrd_3'].iterrows():
        neq_str = '{} != {}'.format(remove_double_quotes(row['VAR1']), remove_double_quotes(row['VAR2']))
        neqls.append(neq_str)
    return neqls

def get_query_heads(ruleH_dfs, ruleHTrue_dfs, hArc_df, colored=True, html=True):
    
    heads = {h: {} for h in set(ruleH_dfs['HEAD'])}
    true_heads = set(ruleHTrue_dfs['HEAD'])
    for _, row in hArc_df.iterrows():
        pos = int(row['POS'])
        h = row['HEAD']
        var = row['VAR']
        heads[h][pos] = var
    heads_htmls = []
    for h in sorted(heads.keys()):
        variables = [remove_double_quotes(heads[h][k]) for k in sorted(heads[h].keys())]
        is_true = h in true_heads
        h_str = '{}{}({})'.format('' if is_true else 'n', head_to_node(h), ','.join(variables))
        h_str = '{}{}{}'.format('' if is_true else 'n', head_to_node(h), '({})'.format(','.join(variables)) if len(variables) > 0 else '')
        if html:
            if colored:
                heads_htmls.append(cstr(h_str, 'green' if is_true else 'red'))
            else:
                heads_htmls.append(cstr(h_str))
        else:
            heads_htmls.append(h_str)
    return heads_htmls

def get_original_query_heads(pw_rel_dfs, colored=True, html=True):
    
    return get_query_heads(pw_rel_dfs['ruleH_1'], pw_rel_dfs['ruleHTrue_1'], pw_rel_dfs['hArc_3'], colored, html)

def get_substituted_query_heads(pw_rel_dfs, colored=True, html=True):
    
    return get_query_heads(pw_rel_dfs['ruleH_1'], pw_rel_dfs['ruleHTrue_1'], pw_rel_dfs['newHArc_3'], colored, html)

def get_query_body_rules(ruleOcc_df, ruleOccTrue_df, arc_df, colored=True, html=True):
    
    rules = {}
    for _, p in ruleOcc_df[['ATOM', 'OCC']].iterrows():
        rules[(p['ATOM'], int(p['OCC']))] = {}
    
    true_rules = set([])
    for _, p in ruleOccTrue_df[['ATOM', 'OCC']].iterrows():
        true_rules.add((p['ATOM'], int(p['OCC'])))
    
    for _, row in arc_df.iterrows():
        rule = (row['ATOM'], int(row['OCC']))
        pos = int(row['POS'])
        var = row['VAR']
        rules[rule][pos] = var
    
    rules_htmls = []
    for (r, occ) in sorted(rules.keys()):
        variables = [remove_double_quotes(rules[(r,occ)][k]) for k in sorted(rules[(r,occ)].keys())]
        is_true = (r,occ) in true_rules
        r_str = '{}{}({})'.format('' if is_true else 'not ', head_to_node(r), ','.join(variables))
        if html:
            if colored:
                rules_htmls.append(cstr(r_str, 'green' if is_true else 'red'))
            else:
                rules_htmls.append(cstr(r_str))
        else:
            rules_htmls.append(r_str)
    
    return rules_htmls

def get_original_query_body_rules(pw_rel_dfs, colored=True, html=True):
    
    return get_query_body_rules(pw_rel_dfs['ruleOcc_2'], pw_rel_dfs['ruleOccTrue_2'], pw_rel_dfs['arc_4'], colored, html)

def get_substituted_query_body_rules(pw_rel_dfs, colored=True, html=True):
    
    return get_query_body_rules(pw_rel_dfs['ruleOcc_2'], pw_rel_dfs['ruleOccTrue_2'], pw_rel_dfs['newArc_4'], colored, html)
    
def print_rewritten_query_string(pw_rel_dfs, html=True, display_q=True, include_ineqls=True):

    heads_htmls = get_substituted_query_heads(pw_rel_dfs, colored=True, html=html)
    rules_htmls = get_substituted_query_body_rules(pw_rel_dfs, colored=True, html=html)
    neqls = neq_rules(pw_rel_dfs)
    if include_ineqls:
        rules_htmls.extend(neqls)
    
    heads_htmls_combined = ' ; '.join(heads_htmls)
    rules_htmls_combined = ', '.join(rules_htmls)
    
    q_str = '{} :- {}.'.format(heads_htmls_combined, rules_htmls_combined)
    
    if html:
        q_str = cstr(q_str)
        if display_q:
            display(html_print(q_str))
    else:
        if display_q:
            print(q_str)
    return q_str
        

def print_explicit_rewritten_query_string(pw_rel_dfs, chain_eq=True, include_ineqls=True, include_eqls=True, display_q=True, html=True):
    
    heads_htmls = get_original_query_heads(pw_rel_dfs, colored=True, html=html)
    rules_htmls = get_original_query_body_rules(pw_rel_dfs, colored=True, html=html)
    eqls = chain_eqs_rules(pw_rel_dfs) if chain_eq else eq_rules(pw_rel_dfs)
    neqls = neq_rules(pw_rel_dfs)
    if include_eqls:
        rules_htmls.extend(eqls)
    if include_ineqls:
        rules_htmls.extend(neqls)
    
    heads_htmls_combined = ' ; '.join(heads_htmls)
    rules_htmls_combined = ', '.join(rules_htmls)
    
    q_str = '{} :- {}.'.format(heads_htmls_combined, rules_htmls_combined)
    if html:
        q_str = cstr(q_str)
        if display_q:
            display(html_print(q_str))
    else:
        if display_q:
            print(q_str)
    return q_str
    
def print_fancy_rewrite(pw_rel_dfs, html=True, display_q=True):
    
    heads_htmls = get_original_query_heads(pw_rel_dfs, colored=True, html=html)
    rules_htmls = get_original_query_body_rules(pw_rel_dfs, colored=True, html=html)
    eq_grps = eq_groups(pw_rel_dfs, single_grps=True)
    eqls = ['='.join(map(remove_double_quotes, grp)) for grp in eq_grps]
    eqls = ['[{}]'.format(eq_grp) for eq_grp in eqls]
    
    heads_htmls_combined = ' ; '.join(heads_htmls)
    rules_htmls_combined = ', '.join(rules_htmls)
    eqls_combined = ''.join(eqls)
    if html:
        eqls_combined = cstr(eqls_combined, 'blue')
    
    q_str = '{} :- {}. % {}'.format(heads_htmls_combined, rules_htmls_combined, eqls_combined)
    
    if html:
        q_str = cstr(q_str)
        if display_q:
            display(html_print(q_str))
    else:
        if display_q:
            print(q_str)
    return q_str
    

def get_query_heads_dict(pw_obj):
    heads = {h: {} for h in set([t[0] for t in pw_obj.rls['ruleH_1']])}
    true_heads = set([t[0] for t in pw_obj.rls['ruleHTrue_1']])
    if 'newHArc_3' in pw_obj.rls:
        for row in pw_obj.rls['newHArc_3']:
            pos, h, var = int(row[1]), row[2], row[0]
            heads[h][pos] = var
    return heads, true_heads
    
def get_query_head_facts(pw_obj, idx=None):
    heads, true_heads = get_query_heads_dict(pw_obj)
    head_strs = []
    for h in sorted(heads.keys()):
        variables = [(heads[h][k]) for k in sorted(heads[h].keys())]
        is_true = h in true_heads
        h_str = '{}{}{}{}.'.format('' if is_true else 'n', head_to_node(h), '' if idx is None else str(idx), '({})'.format(','.join(variables)) if len(variables) > 0 else '')
        head_strs.append(h_str)
    return head_strs

def get_edge_facts(pw, edge_rel_idx: int=None):
    pw_objs = [pw]
    _,_,pw_objs = rel_slicer(None, None, pw_objs, rels_to_use=['e_2', 'ne_2'])
    if edge_rel_idx is not None:
        _, pw_objs, _ = rel_name_remapper(None, pw_objs, None, rel_name_map={'e_2': 'e{}_2'.format(str(edge_rel_idx)),
                                                                             'ne_2': 'ne{}_2'.format(str(edge_rel_idx))})
    return PWEExport.export_as_asp_facts(pw_objs, include_pw_ids=False)

def are_equivalent_patterns(pw1, pw2, eq_check_encoding):
    pw1_edge_facts = get_edge_facts(pw1, 1)
    pw2_edge_facts = get_edge_facts(pw2, 2)
    pw1_head_facts = get_query_head_facts(pw1, idx=1)
    pw2_head_facts = get_query_head_facts(pw2, idx=2)
    #print(pw1_edge_facts, pw2_edge_facts, pw1_head_facts, pw2_head_facts)
    asp_out, _ = run_clingo(eq_check_encoding.split('\n')+pw1_edge_facts+pw2_edge_facts+pw1_head_facts+pw2_head_facts, num_solutions=1)
    #print(asp_out)
    _,_,eq_check_pws = load_worlds(asp_out, silent=True)
    
    return len(eq_check_pws) >= 1


def get_equivalent_sets(objs, match_func):
    
    sets = []
    curr_iter_set = list(range(len(objs)))
    while len(curr_iter_set) > 0:
        next_iter_set = []
        seed = curr_iter_set[0]
        curr_set = {objs[seed]}
        for c in curr_iter_set[1:]:
            if match_func(objs[seed], objs[c]):
                curr_set.add(objs[c])
            else:
                next_iter_set.append(c)
        sets.append(curr_set)
        curr_iter_set = next_iter_set
    return sets

def get_eqs_set(pw_rel_dfs, pw_id):
    eq_df = pw_rel_dfs['eqOrd_3']
    pw_eq_df = eq_df[eq_df['pw'] == pw_id]
    pw_eq_df = pw_eq_df[['VAR1', 'VAR2']]
    return set((r['VAR1'], r['VAR2']) for _,r in pw_eq_df.iterrows())


from itertools import chain, combinations
def powerset(iterable):
    "powerset([1,2,3]) --> () (1,) (2,) (3,) (1,2) (1,3) (2,3) (1,2,3)"
    s = list(iterable)
    return chain.from_iterable(combinations(s, r) for r in range(len(s)+1))