import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Circle

# ================= 配置参数 =================
# 场地参数 (RoboLeague)
FIELD_LENGTH = 14.0
FIELD_WIDTH = 9.0
PENALTY_DIST = 2.1
GOAL_WIDTH = 2.6
CIRCLE_RADIUS = 1.5
PENALTY_AREA_LEN = 3.0
PENALTY_AREA_WID = 6.0
GOAL_AREA_LEN = 1.0
GOAL_AREA_WID = 4.0

# 物理限制
V_X_MAX = 2.0
V_Y_MAX = 1.5
OMEGA_MAX = 1.5

# 绘图分辨率
GRID_RES = 0.2

# ======= 视觉层叠：底色(符号) + 热力图(强弱) =======
POS_COLOR = "#f6b3b3"   # 淡红
NEG_COLOR = "#bcd9ff"   # 淡蓝
CTRL_ALPHA = 0.28       # 底色透明度
HEAT_ALPHA = 0.55       # 热力图透明度（叠在底色上）
ZERO_LINE_W = 2.2

# ======= vx 手柄可点性 =======
VX_VIS_MIN = 0.40       # 显示用最小箭头长度（vx 可为 0，但仍能点到手柄）
TIP_PICK_R = 0.90       # 手柄点击半径
BODY_PICK_R = 0.80
BALL_PICK_R = 0.45

# 点击优先级：越小越优先（dist + bias）
BIAS_TIP  = 0.00
BIAS_BODY = 0.08
BIAS_BALL = 0.25

class RobotSoccerSim:
    def __init__(self):
        self.fig, self.ax = plt.subplots(figsize=(12, 8))
        plt.subplots_adjust(bottom=0.08, top=0.92)

        # 初始化状态: [x, y, theta(rad), vx(m/s)]
        self.team_a = [[-3, 0, 0, 1.0], [-6, 3, 0, 1.0]]
        self.team_b = [[3, 0, np.pi, 1.0], [6, -3, np.pi, 1.0]]
        self.ball = [0, 0]

        self.default_vy = 0.0
        self.default_omega = 0.0

        # 拖拽状态
        self.dragging_obj = None
        self.drag_idx = -1

        # 记录最后鼠标位置（用于键盘操作，避免 event.xdata None）
        self.last_mouse_xy = (0.0, 0.0)

        # 生成网格（略超出球场；区域mask只取球场内）
        x = np.arange(-FIELD_LENGTH/2 - 2, FIELD_LENGTH/2 + 2, GRID_RES)
        y = np.arange(-FIELD_WIDTH/2 - 2, FIELD_WIDTH/2 + 2, GRID_RES)
        self.X, self.Y = np.meshgrid(x, y)

        # 事件绑定
        self.fig.canvas.mpl_connect('button_press_event', self.on_click)
        self.fig.canvas.mpl_connect('button_release_event', self.on_release)
        self.fig.canvas.mpl_connect('motion_notify_event', self.on_motion)
        self.fig.canvas.mpl_connect('key_press_event', self.on_key)

        self.draw_scene()
        plt.show()

    @staticmethod
    def mirror_x(xmin, xmax):
        # 关于 x=0 镜像： [xmin, xmax] -> [-xmax, -xmin]
        return (-xmax, -xmin)


    # ----------------- 足球场绘制 -----------------
    def draw_field(self):
        self.ax.add_patch(Rectangle((-FIELD_LENGTH/2, -FIELD_WIDTH/2),
                                    FIELD_LENGTH, FIELD_WIDTH,
                                    fill=False, color='black', linewidth=2, zorder=6))
        self.ax.plot([0, 0], [-FIELD_WIDTH/2, FIELD_WIDTH/2], color='black', zorder=6)
        self.ax.add_patch(Circle((0, 0), CIRCLE_RADIUS, fill=False, color='black', zorder=6))
        self.ax.add_patch(Circle((0, 0), 0.1, color='black', zorder=6))

        # 左右半场用正宽度矩形，避免 sign*width 负宽度
        for side in [-1, 1]:
            if side == -1:
                pa_x0 = -FIELD_LENGTH/2
                ga_x0 = -FIELD_LENGTH/2
                gx0 = -FIELD_LENGTH/2 - 0.5
            else:
                pa_x0 = FIELD_LENGTH/2 - PENALTY_AREA_LEN
                ga_x0 = FIELD_LENGTH/2 - GOAL_AREA_LEN
                gx0 = FIELD_LENGTH/2

            self.ax.add_patch(Rectangle((pa_x0, -PENALTY_AREA_WID/2),
                                        PENALTY_AREA_LEN, PENALTY_AREA_WID,
                                        fill=False, color='black', zorder=6))
            self.ax.add_patch(Rectangle((ga_x0, -GOAL_AREA_WID/2),
                                        GOAL_AREA_LEN, GOAL_AREA_WID,
                                        fill=False, color='black', zorder=6))
            px = side * (FIELD_LENGTH/2 - PENALTY_DIST)
            self.ax.add_patch(Circle((px, 0), 0.1, color='black', zorder=6))

            self.ax.add_patch(Rectangle((gx0, -GOAL_WIDTH/2),
                                        0.5, GOAL_WIDTH,
                                        color='gray', alpha=0.5, zorder=6))


    # ----------------- 势场计算 -----------------
    def calculate_score_field(self):
        Z_total = np.zeros_like(self.X)

        def get_local_score(rx, ry, r_theta, vx, vy, w, ball_x, ball_y):
            dx_global = self.X - rx - 0.5*vx*np.cos(r_theta) + 0.5*vy*np.sin(r_theta)
            dy_global = self.Y - ry - 0.5*vx*np.sin(r_theta) - 0.5*vy*np.cos(r_theta)

            cos_t = np.cos(-r_theta)
            sin_t = np.sin(-r_theta)

            local_x = dx_global * cos_t - dy_global * sin_t
            local_y = dx_global * sin_t + dy_global * cos_t

            dbx_global = ball_x - rx - 0.5*vx*np.cos(r_theta) + 0.5*vy*np.sin(r_theta)
            dby_global = ball_y - ry - 0.5*vx*np.sin(r_theta) - 0.5*vy*np.cos(r_theta)
            local_a = dbx_global * cos_t - dby_global * sin_t
            local_b = dbx_global * sin_t + dby_global * cos_t

            dist_term = np.sqrt(local_x**2 + local_y**2) / V_X_MAX

            denom_angle = np.sqrt(1e-6 + (local_x**2 + local_y**2) * ((local_a - local_x)**2 + (local_b - local_y)**2))
            angle_input = (local_x*(local_a - local_x) + local_y*(local_b - local_y)) / denom_angle
            angle_input = np.clip(angle_input, -1.0, 1.0)

            term1_input = local_x / np.sqrt(1e-6 + local_x**2 + local_y**2)
            term1_input = np.clip(term1_input, -1.0, 1.0)

            angle_term = (np.arccos(term1_input) + np.arccos(angle_input)) / OMEGA_MAX
            exponent = -(dist_term + angle_term)
            return np.exp(exponent)

        for p in self.team_a:
            Z_total += get_local_score(p[0], p[1], p[2], p[3],
                                       self.default_vy, self.default_omega,
                                       self.ball[0], self.ball[1])

        for p in self.team_b:
            Z_total -= get_local_score(p[0], p[1], p[2], p[3],
                                       self.default_vy, self.default_omega,
                                       self.ball[0], self.ball[1])

        return Z_total

    # ----------------- 绘制场景（含实时统计） -----------------
    def draw_scene(self):
        self.ax.clear()

        Z = self.calculate_score_field()

        # ========== 1) 底色：>0 淡红，<0 淡蓝 ==========
        zmin, zmax = float(np.min(Z)), float(np.max(Z))
        eps = 1e-9
        self.ax.contourf(self.X, self.Y, Z,
                         levels=[zmin - eps, 0, zmax + eps],
                         colors=[NEG_COLOR, POS_COLOR],
                         alpha=CTRL_ALPHA, zorder=0)

        # ========== 2) 强弱热力图：叠加上色 ==========
        max_val = max(abs(zmin), abs(zmax), 0.1)
        levels = np.linspace(-max_val, max_val, 60)
        self.ax.contourf(self.X, self.Y, Z,
                         levels=levels, cmap='RdBu_r',
                         alpha=HEAT_ALPHA, extend='both', zorder=1)

        # 0 分界线
        self.ax.contour(self.X, self.Y, Z, levels=[0],
                        colors='white', linewidths=ZERO_LINE_W, zorder=2)

        # ========== 3) 球场线 ==========
        self.draw_field()

        # ========== 4) 球 & 机器人 ==========
        self.ax.plot(self.ball[0], self.ball[1], 'o',
                     color='gold', markeredgecolor='black',
                     markersize=14, zorder=30)

        for p in self.team_a:
            self.draw_robot(p, 'darkred')
        for p in self.team_b:
            self.draw_robot(p, 'darkblue')

        # 装饰
        self.ax.set_title("Reward Field + 22-zone weighted territorial score\n"
                          "Drag body/tip, A/B add, D delete (tip>body>ball picking)")
        self.ax.set_xlim(-FIELD_LENGTH/2 - 2, FIELD_LENGTH/2 + 2)
        self.ax.set_ylim(-FIELD_WIDTH/2 - 2, FIELD_WIDTH/2 + 2)
        self.ax.set_aspect('equal')
        self.fig.canvas.draw()

    # ----------------- 机器人绘制与可点 tip 位置 -----------------
    def tip_pos_visual(self, state):
        x, y, theta, vx = state
        L = max(float(vx), VX_VIS_MIN)  # 显示用长度
        return x + L*np.cos(theta), y + L*np.sin(theta)

    def draw_robot(self, state, color):
        x, y, theta, vx = state

        self.ax.plot(x, y, 'o', color=color,
                     markeredgecolor='white', markeredgewidth=2,
                     markersize=12, zorder=20)

        tip_x, tip_y = self.tip_pos_visual(state)
        self.ax.plot([x, tip_x], [y, tip_y], color='white', lw=3, zorder=21)
        self.ax.plot(tip_x, tip_y, 'o', color='white',
                     markersize=7.5, markeredgecolor='black', zorder=22)

    # ----------------- 选中逻辑 -----------------
    def find_nearest(self, x, y):
        cands = []

        # tip 优先
        for i, p in enumerate(self.team_a):
            tx, ty = self.tip_pos_visual(p)
            d = np.hypot(tx-x, ty-y)
            if d < TIP_PICK_R:
                cands.append((d + BIAS_TIP, 'A_arrow', i))

        for i, p in enumerate(self.team_b):
            tx, ty = self.tip_pos_visual(p)
            d = np.hypot(tx-x, ty-y)
            if d < TIP_PICK_R:
                cands.append((d + BIAS_TIP, 'B_arrow', i))

        # body 其次
        for i, p in enumerate(self.team_a):
            d = np.hypot(p[0]-x, p[1]-y)
            if d < BODY_PICK_R:
                cands.append((d + BIAS_BODY, 'A_body', i))

        for i, p in enumerate(self.team_b):
            d = np.hypot(p[0]-x, p[1]-y)
            if d < BODY_PICK_R:
                cands.append((d + BIAS_BODY, 'B_body', i))

        # ball 最后
        db = np.hypot(self.ball[0]-x, self.ball[1]-y)
        if db < BALL_PICK_R:
            cands.append((db + BIAS_BALL, 'ball', 0))

        if not cands:
            return None, -1
        cands.sort(key=lambda t: t[0])
        return cands[0][1], cands[0][2]

    # ----------------- 事件 -----------------
    def on_click(self, event):
        if event.inaxes != self.ax or event.xdata is None or event.ydata is None:
            return
        self.last_mouse_xy = (event.xdata, event.ydata)

        t, idx = self.find_nearest(event.xdata, event.ydata)
        if t is None:
            return
        self.dragging_obj = t
        self.drag_idx = idx

    def on_release(self, event):
        self.dragging_obj = None
        self.drag_idx = -1

    def on_motion(self, event):
        if event.inaxes != self.ax or event.xdata is None or event.ydata is None:
            return

        self.last_mouse_xy = (event.xdata, event.ydata)

        if not self.dragging_obj:
            return

        x, y = event.xdata, event.ydata

        if self.dragging_obj == 'ball':
            self.ball = [x, y]

        elif self.dragging_obj == 'A_body':
            self.team_a[self.drag_idx][0] = x
            self.team_a[self.drag_idx][1] = y

        elif self.dragging_obj == 'B_body':
            self.team_b[self.drag_idx][0] = x
            self.team_b[self.drag_idx][1] = y

        elif self.dragging_obj == 'A_arrow':
            rx, ry, theta, vx = self.team_a[self.drag_idx]
            dx, dy = x - rx, y - ry
            raw_dist = float(np.hypot(dx, dy))

            # vx 最大截断到 V_X_MAX，允许 vx=0（可点性由 VX_VIS_MIN 保证）
            self.team_a[self.drag_idx][3] = float(np.clip(raw_dist, 0.0, V_X_MAX))
            if raw_dist > 1e-4:
                self.team_a[self.drag_idx][2] = float(np.arctan2(dy, dx))

        elif self.dragging_obj == 'B_arrow':
            rx, ry, theta, vx = self.team_b[self.drag_idx]
            dx, dy = x - rx, y - ry
            raw_dist = float(np.hypot(dx, dy))

            self.team_b[self.drag_idx][3] = float(np.clip(raw_dist, 0.0, V_X_MAX))
            if raw_dist > 1e-4:
                self.team_b[self.drag_idx][2] = float(np.arctan2(dy, dx))

        self.draw_scene()

    def on_key(self, event):
        x, y = self.last_mouse_xy

        if event.key == 'a':
            self.team_a.append([x, y, 0.0, 1.0])
            self.draw_scene()
            return

        if event.key == 'b':
            self.team_b.append([x, y, np.pi, 1.0])
            self.draw_scene()
            return

        if event.key == 'd':
            t, idx = self.find_nearest(x, y)
            if t is None or t == 'ball':
                return
            if t.startswith('A'):
                self.team_a.pop(idx)
            elif t.startswith('B'):
                self.team_b.pop(idx)
            self.draw_scene()
            return

        return


if __name__ == "__main__":
    app = RobotSoccerSim()