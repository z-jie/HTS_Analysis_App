import sys
import os
import numpy as np
import pandas as pd
import tkinter as tk
from tkinter import ttk, filedialog, messagebox, simpledialog
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.figure import Figure
from matplotlib.backends.backend_agg import FigureCanvasAgg 
from matplotlib import font_manager
import matplotlib.pyplot as plt
from PIL import Image, ImageTk
import traceback 

# =========================================================================
# [核心修复] 显式导入 matplotlib 后端，防止 PyInstaller 打包丢失
# =========================================================================
import matplotlib.backends.backend_pdf
import matplotlib.backends.backend_svg
import matplotlib.backends.backend_ps

# =========================================================================
# 资源路径处理函数
# =========================================================================
def resource_path(relative_path):
    """ 获取资源的绝对路径 """
    try:
        base_path = sys._MEIPASS
    except Exception:
        base_path = os.path.abspath(".")
    return os.path.join(base_path, relative_path)

# =========================================================================
# 自动配置中文字体
# =========================================================================
def configure_chinese_font():
    """ 尝试找到并设置可用的中文字体 """
    system_fonts = font_manager.findSystemFonts(fontpaths=None, fontext='ttf')
    preferred_fonts = ['SimHei', 'Microsoft YaHei', 'SimSun', 'Arial Unicode MS', 'DengXian']
    
    found_font = None
    for font_name in preferred_fonts:
        try:
            plt.rcParams['font.sans-serif'] = [font_name]
            plt.rcParams['axes.unicode_minus'] = False 
            found_font = font_name
            break
        except:
            continue
            
    if not found_font:
        for font_path in system_fonts:
            try:
                font_prop = font_manager.FontProperties(fname=font_path)
                if any(x in font_prop.get_name().lower() for x in ['hei', 'song', 'kai', 'ming']):
                    plt.rcParams['font.sans-serif'] = [font_prop.get_name()]
                    plt.rcParams['axes.unicode_minus'] = False
                    break
            except:
                continue

class HTSAnalysisApp:
    def __init__(self, root):
        configure_chinese_font() 
        
        self.root = root
        self.root.title("高温超导磁悬浮悬浮导向力分析系统 V3.9")
        self.root.geometry("1200x850") 
        
        self.data_store = []
        self.colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', 
                       '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf']
        self.next_color_idx = 0
        self.current_context_ax = None 

        self.setup_ui()

    def setup_ui(self):
        # 1. 顶部 Header
        header_frame = tk.Frame(self.root, bg="white", height=80)
        header_frame.pack(side=tk.TOP, fill=tk.X)
        header_frame.pack_propagate(False)

        try:
            self.photo1 = self.load_keep_ratio_image("西南交通大学.png", 200, 60)
            if self.photo1:
                lbl_logo1 = tk.Label(header_frame, image=self.photo1, bg="white")
                lbl_logo1.pack(side=tk.LEFT, padx=20)

            self.photo2 = self.load_keep_ratio_image("超导磁浮.jpg", 200, 60)
            if self.photo2:
                lbl_logo2 = tk.Label(header_frame, image=self.photo2, bg="white")
                lbl_logo2.pack(side=tk.RIGHT, padx=20)
        except Exception as e:
            print(f"Warning: Logo loading failed. {e}")

        lbl_title = tk.Label(header_frame, text="高温超导磁悬浮悬浮导向力分析系统", 
                             font=("黑体", 20, "bold"), fg="#0072BD", bg="white")
        lbl_title.pack(side=tk.TOP, pady=20)

        # 2. 主内容区域
        main_content = ttk.PanedWindow(self.root, orient=tk.HORIZONTAL)
        main_content.pack(fill=tk.BOTH, expand=True)

        # --- 左侧控制面板 ---
        left_panel = ttk.Frame(main_content, width=300)
        main_content.add(left_panel, weight=0) 

        # 参数输入区
        param_group = ttk.LabelFrame(left_panel, text="控制面板")
        param_group.pack(fill=tk.X, padx=10, pady=10)

        self.params = {}
        defaults = [("场冷高度 (mm):", "30", "FCH"), 
                    ("工作高度 (mm):", "10", "WH"),
                    ("横向偏移 (mm):", "15", "OFFSET"),
                    ("块材长度 (m):", "0.032", "Length")]
        
        for i, (label_text, default_val, key) in enumerate(defaults):
            lbl = ttk.Label(param_group, text=label_text)
            lbl.grid(row=i, column=0, padx=5, pady=5, sticky="e")
            entry = ttk.Entry(param_group, width=10)
            entry.insert(0, default_val)
            entry.grid(row=i, column=1, padx=5, pady=5)
            self.params[key] = entry

        # 按钮区
        btn_frame = ttk.Frame(left_panel)
        btn_frame.pack(fill=tk.X, padx=10, pady=5)

        ttk.Button(btn_frame, text="导入仿真数据 (.xlsx/.csv)", command=self.import_sim_data).pack(fill=tk.X, pady=2)
        ttk.Button(btn_frame, text="导入实验数据 (.xlsx/.csv)", command=self.import_exp_data).pack(fill=tk.X, pady=2)
        
        ttk.Separator(left_panel, orient='horizontal').pack(fill=tk.X, padx=10, pady=10)

        # 工况列表
        list_frame = ttk.LabelFrame(left_panel, text="工况图层管理")
        list_frame.pack(fill=tk.BOTH, expand=True, padx=10, pady=5)
        
        self.canvas_list = tk.Canvas(list_frame, bg="white")
        self.scrollbar = ttk.Scrollbar(list_frame, orient="vertical", command=self.canvas_list.yview)
        self.scrollable_frame = ttk.Frame(self.canvas_list)

        self.scrollable_frame.bind(
            "<Configure>",
            lambda e: self.canvas_list.configure(scrollregion=self.canvas_list.bbox("all"))
        )
        self.canvas_list.create_window((0, 0), window=self.scrollable_frame, anchor="nw")
        self.canvas_list.configure(yscrollcommand=self.scrollbar.set)

        self.canvas_list.pack(side="left", fill="both", expand=True)
        self.scrollbar.pack(side="right", fill="y")

        # 底部按钮
        bottom_frame = ttk.Frame(left_panel)
        bottom_frame.pack(fill=tk.X, padx=10, pady=10, side=tk.BOTTOM)
        
        ttk.Button(bottom_frame, text="删除选中工况", command=self.delete_layer).pack(fill=tk.X, pady=2)
        ttk.Button(bottom_frame, text="导出处理后数据", command=self.export_data).pack(fill=tk.X, pady=2)
        
        h_btn_frame = ttk.Frame(bottom_frame)
        h_btn_frame.pack(fill=tk.X, pady=5)
        ttk.Button(h_btn_frame, text="帮助/关于", width=12, command=self.show_help).pack(side=tk.LEFT)
        ttk.Button(h_btn_frame, text="清空所有", width=12, command=self.clear_all).pack(side=tk.RIGHT)

        # --- 右侧绘图区域 ---
        right_panel = ttk.Frame(main_content)
        main_content.add(right_panel, weight=4)

        self.fig = Figure(figsize=(8, 6), dpi=100)
        self.fig.subplots_adjust(left=0.1, bottom=0.1, right=0.95, top=0.9, wspace=0.3, hspace=0.4)
        
        self.ax1 = self.fig.add_subplot(221)
        self.ax2 = self.fig.add_subplot(222)
        self.ax3 = self.fig.add_subplot(223)
        self.ax4 = self.fig.add_subplot(224)
        
        self.setup_axes()

        self.canvas_plot = FigureCanvasTkAgg(self.fig, master=right_panel)
        self.canvas_plot.draw()
        
        toolbar = NavigationToolbar2Tk(self.canvas_plot, right_panel)
        toolbar.update()
        self.canvas_plot.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=True)

        self.canvas_plot.mpl_connect('scroll_event', self.on_scroll)
        self.canvas_plot.mpl_connect('button_press_event', self.on_mouse_click)

    def load_keep_ratio_image(self, filename, max_w, max_h):
        path = resource_path(filename)
        if not os.path.exists(path): return None
        img = Image.open(path)
        img.thumbnail((max_w, max_h), Image.Resampling.LANCZOS) 
        return ImageTk.PhotoImage(img)

    def on_scroll(self, event):
        ax = event.inaxes
        if ax is None: return
        
        base_scale = 1.1
        if event.button == 'up': scale_factor = 1 / base_scale
        elif event.button == 'down': scale_factor = base_scale
        else: return

        x_min, x_max = ax.get_xlim()
        y_min, y_max = ax.get_ylim()
        
        x_data = event.xdata
        y_data = event.ydata
        if x_data is None or y_data is None: return

        new_width = (x_max - x_min) * scale_factor
        new_height = (y_max - y_min) * scale_factor
        
        rel_x = (x_max - x_data) / (x_max - x_min)
        rel_y = (y_max - y_data) / (y_max - y_min)
        
        ax.set_xlim([x_data - new_width * (1 - rel_x), x_data + new_width * rel_x])
        ax.set_ylim([y_data - new_height * (1 - rel_y), y_data + new_height * rel_y])
        self.canvas_plot.draw()

    def on_mouse_click(self, event):
        if event.button == 3 and event.inaxes: 
            self.current_context_ax = event.inaxes
            menu = tk.Menu(self.root, tearoff=0)
            menu.add_command(label="另存为图片 (Save Plot)", command=self.save_current_ax)
            menu.tk_popup(event.guiEvent.x_root, event.guiEvent.y_root)

    def save_current_ax(self):
        ax = self.current_context_ax
        if not ax: return
        
        filepath = filedialog.asksaveasfilename(defaultextension=".png", 
                                                filetypes=[("PNG Image", "*.png"), ("JPEG Image", "*.jpg"), ("PDF Document", "*.pdf")])
        if not filepath: return

        try:
            temp_fig = Figure(figsize=(6, 5), dpi=150)
            temp_ax = temp_fig.add_subplot(111)
            
            temp_ax.set_title(ax.get_title())
            temp_ax.set_xlabel(ax.get_xlabel())
            temp_ax.set_ylabel(ax.get_ylabel())
            temp_ax.set_xlim(ax.get_xlim())
            temp_ax.set_ylim(ax.get_ylim())
            temp_ax.grid(True, linestyle='--', alpha=0.6)
            
            for line in ax.get_lines():
                temp_ax.plot(line.get_xdata(), line.get_ydata(), 
                             label=line.get_label(), 
                             color=line.get_color(),
                             linestyle=line.get_linestyle(),
                             linewidth=line.get_linewidth(),
                             marker=line.get_marker())
            
            if ax.get_legend():
                temp_ax.legend()
            
            canvas = FigureCanvasAgg(temp_fig)
            temp_fig.savefig(filepath)
            messagebox.showinfo("成功", "图片保存成功！")
        except Exception as e:
            messagebox.showerror("保存失败", f"保存图片时出错:\n{str(e)}")

    def setup_axes(self):
        ax_list = [
            (self.ax1, '仿真-原始悬浮力数据', '累积位移 / 时间步', '悬浮力 (N)'),
            (self.ax2, '仿真-原始导向力数据', '累积位移 / 时间步', '导向力 (N)'),
            (self.ax3, '悬浮特性 (磁滞回线)', '悬浮间隙 Gap (mm)', '悬浮力 (N)'),
            (self.ax4, '导向特性 (磁滞回线)', '横向位移 (mm)', '导向力 (N)')
        ]
        for ax, title, xl, yl in ax_list:
            ax.set_title(title, fontsize=10, fontweight='bold')
            ax.set_xlabel(xl, fontsize=9)
            ax.set_ylabel(yl, fontsize=9)
            ax.grid(True, linestyle='--', alpha=0.6)

    def safe_read_data(self, filepath, min_cols):
        try:
            df = None
            if filepath.endswith(('.xlsx', '.xls')):
                try:
                    df = pd.read_excel(filepath, header=None, usecols="A:H")
                except Exception:
                    pass
            
            if df is None:
                if filepath.endswith(('.xlsx', '.xls')):
                    df = pd.read_excel(filepath, header=None)
                else:
                    df = pd.read_csv(filepath, header=None)
            
            df = df.apply(pd.to_numeric, errors='coerce')
            df.dropna(how='all', inplace=True)
            
            df.columns = range(df.shape[1])
            
            current_cols = df.shape[1]
            if current_cols < min_cols:
                for i in range(current_cols, min_cols):
                    df[i] = np.nan
            
            df = df.sort_index(axis=1)
            df = df.iloc[:, :min_cols]
            
            return df
        except Exception as e:
            messagebox.showerror("读取错误", f"无法读取文件: {str(e)}")
            return None

    def is_valid_column(self, series):
        if series.isna().all(): return False
        valid_data = series.dropna()
        if (valid_data.abs() < 1e-9).all(): return False
        return True

    def process_sim_data(self, df, p):
        X_sim = df.iloc[:, 0].values
        
        if self.is_valid_column(df.iloc[:, 1]):
            F_lev = df.iloc[:, 1].values * p['Length']
            has_lev = True
        else:
            F_lev = None
            has_lev = False
            
        if self.is_valid_column(df.iloc[:, 2]):
            F_gui = df.iloc[:, 2].values * p['Length']
            has_gui = True
        else:
            F_gui = None
            has_gui = False
            
        H_drop = p['FCH'] - p['WH']
        T_start = 4 * H_drop
        
        res = {'Gap': None, 'Fz': None, 'Lat': None, 'Fy': None, 
               'RawX': X_sim, 'RawFz': F_lev, 'RawFy': F_gui,
               'H_drop': H_drop, 'T_start': T_start, 
               'MaxGap': p['FCH'], 'MaxOffset': p['OFFSET']}
        
        if has_lev:
            idx_des = (X_sim >= 0) & (X_sim <= H_drop)
            idx_asc = (X_sim > H_drop) & (X_sim <= 2 * H_drop)
            
            gap_des = p['FCH'] - X_sim[idx_des]
            gap_asc = p['WH'] + (X_sim[idx_asc] - H_drop)
            
            if np.any(gap_des < -1.0): 
                 messagebox.showwarning("参数警告", f"检测到悬浮间隙为负值。请检查场冷高度设置。")

            res['Gap'] = np.concatenate([gap_des, gap_asc])
            res['Fz'] = np.concatenate([F_lev[idx_des], F_lev[idx_asc]])
            
        if has_gui:
            range_r = [T_start, T_start + p['OFFSET']]
            range_l = [T_start + p['OFFSET'], T_start + 3*p['OFFSET']]
            range_b = [T_start + 3*p['OFFSET'], T_start + 4*p['OFFSET']]
            
            idx_r = (X_sim >= range_r[0]) & (X_sim <= range_r[1])
            idx_l = (X_sim > range_l[0]) & (X_sim <= range_l[1])
            idx_b = (X_sim > range_b[0]) & (X_sim <= range_b[1])
            
            lat_r = X_sim[idx_r] - T_start
            lat_l = p['OFFSET'] - (X_sim[idx_l] - (T_start + p['OFFSET']))
            lat_b = -p['OFFSET'] + (X_sim[idx_b] - (T_start + 3*p['OFFSET']))
            
            res['Lat'] = np.concatenate([lat_r, lat_l, lat_b])
            res['Fy'] = np.concatenate([F_gui[idx_r], F_gui[idx_l], F_gui[idx_b]])
            
        return res

    def get_ui_params(self):
        try:
            return {
                'FCH': float(self.params['FCH'].get()),
                'WH': float(self.params['WH'].get()),
                'OFFSET': float(self.params['OFFSET'].get()),
                'Length': float(self.params['Length'].get())
            }
        except ValueError:
            messagebox.showerror("错误", "请输入有效的数字参数")
            return None

    def add_to_list(self, name, type_str, proc_data, raw_data, params_data):
        color = self.colors[self.next_color_idx % len(self.colors)]
        self.next_color_idx += 1
        var = tk.BooleanVar(value=True)
        var.trace_add("write", lambda *args: self.update_plots())
        
        item = {
            'name': name, 'type': type_str, 'proc': proc_data, 'raw': raw_data,
            'params': params_data, 'color': color, 'var': var, 'widget_ref': None 
        }
        self.data_store.append(item)
        
        frame = ttk.Frame(self.scrollable_frame)
        frame.pack(fill=tk.X, pady=2)
        chk = tk.Checkbutton(frame, variable=var, bg="white")
        chk.pack(side=tk.LEFT)
        lbl_color = tk.Label(frame, bg=color, width=2)
        lbl_color.pack(side=tk.LEFT, padx=5)
        tk.Label(frame, text=f"[{type_str}] {name}").pack(side=tk.LEFT)
        item['widget_ref'] = frame
        self.update_plots()

    def import_sim_data(self):
        try:
            p = self.get_ui_params()
            if not p: return
            filepath = filedialog.askopenfilename(filetypes=[("Excel Files", "*.xlsx;*.xls"), ("CSV Files", "*.csv")])
            if not filepath: return
            default_name = f"Sim_FC{int(p['FCH'])}"
            name = simpledialog.askstring("工况命名", "请输入该仿真工况名称:", initialvalue=default_name)
            if not name: return
            
            df = self.safe_read_data(filepath, 3)
            if df is None: return
            if df.shape[0] < 2:
                 messagebox.showerror("数据错误", "文件有效数据行数不足，请检查文件格式。")
                 return

            proc_data = self.process_sim_data(df, p)
            raw_data = {'X': proc_data['RawX'], 'Fz': proc_data['RawFz'], 'Fy': proc_data['RawFy']}
            params_data = {'H_drop': proc_data['H_drop'], 'T_start': proc_data['T_start']}
            self.add_to_list(name, 'Sim', proc_data, raw_data, params_data)
        except Exception as e:
            messagebox.showerror("导入失败", f"仿真数据处理出错:\n{str(e)}\n{traceback.format_exc()}")

    def import_exp_data(self):
        try:
            filepath = filedialog.askopenfilename(filetypes=[("Excel Files", "*.xlsx;*.xls"), ("CSV Files", "*.csv")])
            if not filepath: return
            name = simpledialog.askstring("工况命名", "请输入该实验工况名称:", initialvalue="Exp_Data")
            if not name: return
            
            df = self.safe_read_data(filepath, 4)
            if df is None: return
            if df.shape[0] < 2:
                 messagebox.showerror("数据错误", "文件有效数据行数不足。")
                 return

            res = {'Gap': None, 'Fz': None, 'Lat': None, 'Fy': None, 'MaxGap':0, 'MaxOffset':0}
            
            if self.is_valid_column(df.iloc[:,0]) and self.is_valid_column(df.iloc[:,1]):
                valid_rows = df.iloc[:,0].notna() & df.iloc[:,1].notna()
                res['Gap'] = df.iloc[:, 0][valid_rows].values
                res['Fz'] = df.iloc[:, 1][valid_rows].values
                if len(res['Gap']) > 0: res['MaxGap'] = np.max(res['Gap'])
                
            if self.is_valid_column(df.iloc[:,2]) and self.is_valid_column(df.iloc[:,3]):
                valid_rows = df.iloc[:,2].notna() & df.iloc[:,3].notna()
                res['Lat'] = df.iloc[:, 2][valid_rows].values
                res['Fy'] = df.iloc[:, 3][valid_rows].values
                if len(res['Lat']) > 0: res['MaxOffset'] = np.max(np.abs(res['Lat']))
            
            self.add_to_list(name, 'Exp', res, None, None)
        except Exception as e:
            messagebox.showerror("导入失败", f"实验数据处理出错:\n{str(e)}\n{traceback.format_exc()}")

    def update_plots(self):
        self.ax1.clear(); self.ax2.clear(); self.ax3.clear(); self.ax4.clear()
        self.setup_axes()
        
        has_sim_lev, has_sim_gui = False, False
        has_proc_lev, has_proc_gui = False, False
        global_max_gap, global_max_offset = 0, 0
        
        for item in self.data_store:
            if not item['var'].get(): continue
            c = item['color']; nm = item['name']
            proc = item['proc']
            
            if proc['Gap'] is not None and len(proc['Gap']) > 0:
                val = np.max(proc['Gap'])
                if 'MaxGap' in proc: val = max(val, proc['MaxGap'])
                global_max_gap = max(global_max_gap, val)
            if proc['Lat'] is not None and len(proc['Lat']) > 0:
                val = np.max(np.abs(proc['Lat']))
                if 'MaxOffset' in proc: val = max(val, proc['MaxOffset'])
                global_max_offset = max(global_max_offset, val)
                
            if item['type'] == 'Sim':
                raw = item['raw']; p = item['params']
                if raw['Fz'] is not None:
                    self.ax1.plot(raw['X'], raw['Fz'], '-', color=c, label=nm)
                    self.ax1.axvline(x=p['H_drop'], color='gray', linestyle=':')
                    has_sim_lev = True
                if raw['Fy'] is not None:
                    self.ax2.plot(raw['X'], raw['Fy'], '-', color=c, label=nm)
                    self.ax2.axvline(x=p['T_start'], color='gray', linestyle=':')
                    has_sim_gui = True
                    
            ls = '--' if item['type'] == 'Exp' else '-'
            if proc['Gap'] is not None and proc['Fz'] is not None:
                self.ax3.plot(proc['Gap'], proc['Fz'], linestyle=ls, color=c, label=nm)
                has_proc_lev = True
            if proc['Lat'] is not None and proc['Fy'] is not None:
                self.ax4.plot(proc['Lat'], proc['Fy'], linestyle=ls, color=c, label=nm)
                has_proc_gui = True
                
        if has_sim_lev: self.ax1.legend()
        if has_sim_gui: self.ax2.legend()
        if has_proc_lev: 
            self.ax3.legend()
            if global_max_gap > 0: self.ax3.set_xlim(0, global_max_gap + 10)
        if has_proc_gui: 
            self.ax4.legend()
            if global_max_offset > 0: 
                axis_margin = 5
                self.ax4.set_xlim(-(global_max_offset + axis_margin), (global_max_offset + axis_margin))
        self.canvas_plot.draw()

    def delete_layer(self):
        if not self.data_store: return
        names = [item['name'] for item in self.data_store]
        del_win = tk.Toplevel(self.root)
        del_win.title("选择要删除的工况")
        del_win.geometry("300x400")
        lb = tk.Listbox(del_win, selectmode=tk.MULTIPLE)
        lb.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)
        for n in names: lb.insert(tk.END, n)
        def confirm_del():
            selections = lb.curselection()
            if not selections: return
            for idx in sorted(selections, reverse=True):
                self.data_store[idx]['widget_ref'].destroy()
                del self.data_store[idx]
            self.update_plots()
            del_win.destroy()
        ttk.Button(del_win, text="确定删除", command=confirm_del).pack(pady=10)

    def clear_all(self):
        if not messagebox.askokcancel("确认", "确定要清空所有数据吗？"): return
        for item in self.data_store: item['widget_ref'].destroy()
        self.data_store = []
        self.next_color_idx = 0
        self.update_plots()

    def export_data(self):
        if not self.data_store:
            messagebox.showinfo("提示", "没有数据可导出")
            return
        export_items = [item for item in self.data_store if item['type'] == 'Sim' and item['var'].get()]
        if not export_items:
            messagebox.showinfo("提示", "请先在列表中勾选需要导出的仿真数据")
            return
        filepath = filedialog.asksaveasfilename(defaultextension=".xlsx", filetypes=[("Excel Files", "*.xlsx")])
        if not filepath: return
        try:
            max_rows = 0
            for item in export_items:
                l1 = len(item['proc']['Gap']) if item['proc']['Gap'] is not None else 0
                l2 = len(item['proc']['Lat']) if item['proc']['Lat'] is not None else 0
                max_rows = max(max_rows, l1, l2)
            if max_rows == 0: max_rows = 1
            data_dict = {}
            for item in export_items:
                nm = item['name']
                gap = item['proc']['Gap'] if item['proc']['Gap'] is not None else []
                fz = item['proc']['Fz'] if item['proc']['Fz'] is not None else []
                lat = item['proc']['Lat'] if item['proc']['Lat'] is not None else []
                fy = item['proc']['Fy'] if item['proc']['Fy'] is not None else []
                col_gap = np.full(max_rows, np.nan); col_gap[:len(gap)] = gap
                col_fz = np.full(max_rows, np.nan); col_fz[:len(fz)] = fz
                col_lat = np.full(max_rows, np.nan); col_lat[:len(lat)] = lat
                col_fy = np.full(max_rows, np.nan); col_fy[:len(fy)] = fy
                data_dict[f'Gap_{nm}'] = col_gap
                data_dict[f'Fz_{nm}'] = col_fz
                data_dict[f'Lat_{nm}'] = col_lat
                data_dict[f'Fy_{nm}'] = col_fy
            df_out = pd.DataFrame(data_dict)
            if os.path.exists(filepath):
                try: os.remove(filepath)
                except: pass
            df_out.to_excel(filepath, index=False)
            messagebox.showinfo("成功", "数据导出成功！")
        except Exception as e:
            messagebox.showerror("导出失败", str(e))

    def show_help(self):
        msg = ("【软件名称】高温超导磁悬浮悬浮导向力分析系统 V3.9\n\n"
               "【开发单位】西南交通大学 高温超导磁悬浮课题组\n"
               "【联系人】1300699158@qq.com\n\n"
               "【软件用途】\n"
               "1. 数据后处理：将 COMSOL 等仿真软件输出的单调时间步长数据，自动转换为符合物理特性的磁滞回线（悬浮力-间隙，导向力-位移）。\n"
               "2. 实验验证：直观对比仿真结果与实验测试数据，辅助校准模型参数。\n"
               "3. 数据管理：支持多工况叠加显示、图层管理及处理后数据导出。\n\n"
               "【使用说明】\n"
               "1. 参数设置：在左侧面板填写场冷高度(FCH)、工作高度(WH)、偏移量(OFFSET)及块材长度。\n"
               "2. 仿真导入：支持 Excel/CSV。软件会自动扫描列数据，若无导向力则只画悬浮曲线。\n"
               "3. 实验导入：支持 Excel/CSV。软件会自动识别 1-2 列和 3-4 列的有效性。\n"
               "4. 导出功能：仅导出在表格中被“勾选”的仿真数据。\n"
               "5. 单图保存：在右侧任意绘图区域点击鼠标右键，选择“另存为图片”。")
        messagebox.showinfo("帮助/关于", msg)

if __name__ == "__main__":
    root = tk.Tk()
    app = HTSAnalysisApp(root)
    root.mainloop()