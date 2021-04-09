import numpy as np
import common
import meshtools as mt
import sys
import timecell

def define_neumann_edge(flag_neumann,extra_neumann,coord,topol,all_edges):

    if (flag_neumann == 'example1_compression_inside'):
        neumann_nodes = []
        neumann_values = []
        neumann_edges = 0

        all_edges = np.asarray(all_edges)

        xmin = 0.35
        xmax = 0.65
        yneu = 2.50

        x = coord[:,0]
        y = coord[:,1]

        for iedge,edge in enumerate(all_edges):
            x1 = x[edge[0]]
            x2 = x[edge[1]]
            y1 = y[edge[0]]
            y2 = y[edge[1]]
            if ((abs(y1-yneu)<=1e-12)and(abs(y2-yneu)<=1e-12)):
                if (((x1>xmin-1e-12)and(x1<xmax+1e-12))and
                    ((x2>xmin-1e-12)and(x2<xmax+1e-12))):
                    neumann_nodes.append(edge+1)
                    neumann_edges = neumann_edges + 1
                    neumann_values.append([0.0,-1.0])

        neumann_nodes = np.asarray(neumann_nodes)
        neumann_values = np.asarray(neumann_values)

    elif (flag_neumann == 'example1_traction_inside'):
        neumann_nodes = []
        neumann_values = []
        neumann_edges = 0

        all_edges = np.asarray(all_edges)

        xmin = 0.35
        xmax = 0.65
        yneu = 2.50

        x = coord[:,0]
        y = coord[:,1]

        for iedge,edge in enumerate(all_edges):
            x1 = x[edge[0]]
            x2 = x[edge[1]]
            y1 = y[edge[0]]
            y2 = y[edge[1]]
            if ((abs(y1-yneu)<=1e-12)and(abs(y2-yneu)<=1e-12)):
                if (((x1>xmin-1e-12)and(x1<xmax+1e-12))and
                    ((x2>xmin-1e-12)and(x2<xmax+1e-12))):
                    neumann_nodes.append(edge+1)
                    neumann_edges = neumann_edges + 1
                    neumann_values.append([0.0,+1.0])

        neumann_nodes = np.asarray(neumann_nodes)
        neumann_values = np.asarray(neumann_values)

    elif (flag_neumann == 'example1_compression_semi_inside'):
        neumann_nodes = []
        neumann_values = []
        neumann_edges = 0

        all_edges = np.asarray(all_edges)

        xmin = 0.35
        xmax = 0.65
        yneu = 2.00

        x = coord[:,0]
        y = coord[:,1]

        for iedge,edge in enumerate(all_edges):
            x1 = x[edge[0]]
            x2 = x[edge[1]]
            y1 = y[edge[0]]
            y2 = y[edge[1]]
            if ((abs(y1-yneu)<=1e-12)and(abs(y2-yneu)<=1e-12)):
                if (((x1>xmin-1e-12)and(x1<xmax+1e-12))and
                    ((x2>xmin-1e-12)and(x2<xmax+1e-12))):
                    neumann_nodes.append(edge+1)
                    neumann_edges = neumann_edges + 1
                    neumann_values.append([0.0,-1.0])

        neumann_nodes = np.asarray(neumann_nodes)
        neumann_values = np.asarray(neumann_values)

    elif (flag_neumann == 'example1_compression_semi_inside_steel'):
        neumann_nodes = []
        neumann_values = []
        neumann_edges = 0

        all_edges = np.asarray(all_edges)

        xmin = 0.35
        xmax = 0.65
        yneu = 2.00

        x = coord[:,0]
        y = coord[:,1]

        for iedge,edge in enumerate(all_edges):
            x1 = x[edge[0]]
            x2 = x[edge[1]]
            y1 = y[edge[0]]
            y2 = y[edge[1]]
            if ((abs(y1-yneu)<=1e-12)and(abs(y2-yneu)<=1e-12)):
                if (((x1>xmin-1e-12)and(x1<xmax+1e-12))and
                    ((x2>xmin-1e-12)and(x2<xmax+1e-12))):
                    neumann_nodes.append(edge+1)
                    neumann_edges = neumann_edges + 1
                    neumann_values.append([0.0,-1e6])

        neumann_nodes = np.asarray(neumann_nodes)
        neumann_values = np.asarray(neumann_values)

    elif (flag_neumann == 'example1_traction_semi_inside'):
        neumann_nodes = []
        neumann_values = []
        neumann_edges = 0

        all_edges = np.asarray(all_edges)

        xmin = 0.35
        xmax = 0.65
        yneu = 2.00

        x = coord[:,0]
        y = coord[:,1]

        for iedge,edge in enumerate(all_edges):
            x1 = x[edge[0]]
            x2 = x[edge[1]]
            y1 = y[edge[0]]
            y2 = y[edge[1]]
            if ((abs(y1-yneu)<=1e-12)and(abs(y2-yneu)<=1e-12)):
                if (((x1>xmin-1e-12)and(x1<xmax+1e-12))and
                    ((x2>xmin-1e-12)and(x2<xmax+1e-12))):
                    neumann_nodes.append(edge+1)
                    neumann_edges = neumann_edges + 1
                    neumann_values.append([0.0,+1.0])

        neumann_nodes = np.asarray(neumann_nodes)
        neumann_values = np.asarray(neumann_values)

    elif (flag_neumann == 'example1_compression'):
        neumann_nodes = []
        neumann_values = []
        neumann_edges = 0

        all_edges = np.asarray(all_edges)

        xmin = 0.0
        xmax = 0.3
        yneu = 2.0

        x = coord[:,0]
        y = coord[:,1]

        for iedge,edge in enumerate(all_edges):
            x1 = x[edge[0]]
            x2 = x[edge[1]]
            y1 = y[edge[0]]
            y2 = y[edge[1]]
            if ((abs(y1-yneu)<=1e-12)and(abs(y2-yneu)<=1e-12)):
                if (((x1>xmin-1e-12)and(x1<xmax+1e-12))and
                    ((x2>xmin-1e-12)and(x2<xmax+1e-12))):
                    neumann_nodes.append(edge+1)
                    neumann_edges = neumann_edges + 1
                    neumann_values.append([0.0,-1.0])

        neumann_nodes = np.asarray(neumann_nodes)
        neumann_values = np.asarray(neumann_values)

    elif (flag_neumann == 'example1_traction'):
        neumann_nodes = []
        neumann_values = []
        neumann_edges = 0

        all_edges = np.asarray(all_edges)

        xmin = 0.0
        xmax = 0.3
        yneu = 2.0

        x = coord[:,0]
        y = coord[:,1]

        for iedge,edge in enumerate(all_edges):
            x1 = x[edge[0]]
            x2 = x[edge[1]]
            y1 = y[edge[0]]
            y2 = y[edge[1]]
            if ((abs(y1-yneu)<=1e-12)and(abs(y2-yneu)<=1e-12)):
                if (((x1>xmin-1e-12)and(x1<xmax+1e-12))and
                    ((x2>xmin-1e-12)and(x2<xmax+1e-12))):
                    neumann_nodes.append(edge+1)
                    neumann_edges = neumann_edges + 1
                    neumann_values.append([0.0,+1.0])

        neumann_nodes = np.asarray(neumann_nodes)
        neumann_values = np.asarray(neumann_values)

    elif (flag_neumann == 'example2'):
        neumann_nodes = []
        neumann_values = []
        neumann_edges = 0

        all_edges = np.asarray(all_edges)

        ymin = 0.0
        ymax = 0.3
        xneu = 2.0

        x = coord[:,0]
        y = coord[:,1]

        for iedge,edge in enumerate(all_edges):
            x1 = x[edge[0]]
            x2 = x[edge[1]]
            y1 = y[edge[0]]
            y2 = y[edge[1]]
            if ((abs(x1-xneu)<=1e-12)and(abs(x2-xneu)<=1e-12)):
                if (((y1>ymin-1e-12)and(y1<ymax+1e-12))and
                    ((y2>ymin-1e-12)and(y2<ymax+1e-12))):
                    neumann_nodes.append(edge+1)
                    neumann_edges = neumann_edges + 1
                    neumann_values.append([0.0,-1.0])

        neumann_nodes = np.asarray(neumann_nodes)
        neumann_values = np.asarray(neumann_values)

    elif (flag_neumann == 'example2_inside'):
        neumann_nodes = []
        neumann_values = []
        neumann_edges = 0

        all_edges = np.asarray(all_edges)

        ymin = 0.35
        ymax = 0.65
        xneu = 2.50

        x = coord[:,0]
        y = coord[:,1]

        for iedge,edge in enumerate(all_edges):
            x1 = x[edge[0]]
            x2 = x[edge[1]]
            y1 = y[edge[0]]
            y2 = y[edge[1]]
            if ((abs(x1-xneu)<=1e-12)and(abs(x2-xneu)<=1e-12)):
                if (((y1>ymin-1e-12)and(y1<ymax+1e-12))and
                    ((y2>ymin-1e-12)and(y2<ymax+1e-12))):
                    neumann_nodes.append(edge+1)
                    neumann_edges = neumann_edges + 1
                    neumann_values.append([0.0,-1.0])

        neumann_nodes = np.asarray(neumann_nodes)
        neumann_values = np.asarray(neumann_values)

    elif (flag_neumann == 'example2_semi_inside'):
        neumann_nodes = []
        neumann_values = []
        neumann_edges = 0

        all_edges = np.asarray(all_edges)

        ymin = 0.6
        ymax = 0.9
        xneu = 2.0

        x = coord[:,0]
        y = coord[:,1]

        for iedge,edge in enumerate(all_edges):
            x1 = x[edge[0]]
            x2 = x[edge[1]]
            y1 = y[edge[0]]
            y2 = y[edge[1]]
            if ((abs(x1-xneu)<=1e-12)and(abs(x2-xneu)<=1e-12)):
                if (((y1>ymin-1e-12)and(y1<ymax+1e-12))and
                    ((y2>ymin-1e-12)and(y2<ymax+1e-12))):
                    neumann_nodes.append(edge+1)
                    neumann_edges = neumann_edges + 1
                    neumann_values.append([0.0,-1.0])

        neumann_nodes = np.asarray(neumann_nodes)
        neumann_values = np.asarray(neumann_values)

    elif (flag_neumann == 'example2_semi_inside_steel'):
        neumann_nodes = []
        neumann_values = []
        neumann_edges = 0

        all_edges = np.asarray(all_edges)

        ymin = 0.6
        ymax = 0.9
        xneu = 2.0

        x = coord[:,0]
        y = coord[:,1]

        for iedge,edge in enumerate(all_edges):
            x1 = x[edge[0]]
            x2 = x[edge[1]]
            y1 = y[edge[0]]
            y2 = y[edge[1]]
            if ((abs(x1-xneu)<=1e-12)and(abs(x2-xneu)<=1e-12)):
                if (((y1>ymin-1e-12)and(y1<ymax+1e-12))and
                    ((y2>ymin-1e-12)and(y2<ymax+1e-12))):
                    neumann_nodes.append(edge+1)
                    neumann_edges = neumann_edges + 1
                    neumann_values.append([0.0,-1e6])

        neumann_nodes = np.asarray(neumann_nodes)
        neumann_values = np.asarray(neumann_values)

    else:
        print( 'Flag ', flag_neumann, ' not supported')

    return neumann_edges, neumann_nodes, neumann_values;
