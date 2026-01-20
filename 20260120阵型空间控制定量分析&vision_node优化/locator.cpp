#include "locator.h"
#include "utils/math.h"
#include "utils/misc.h"
#include "brain.h"
#include "utils/print.h"
#include "brain_tree.h"

#include <algorithm>

#define REGISTER_LOCATOR_BUILDER(Name)     \
    factory.registerBuilder<Name>(         \
        #Name,                             \
        [brain](const string &name, const NodeConfig &config) { return make_unique<Name>(name, config, brain); });

void RegisterLocatorNodes(BT::BehaviorTreeFactory &factory, Brain* brain)
{
    REGISTER_LOCATOR_BUILDER(SelfLocate);
    REGISTER_LOCATOR_BUILDER(SelfLocateEnterField);
    REGISTER_LOCATOR_BUILDER(SelfLocate1M);
    REGISTER_LOCATOR_BUILDER(SelfLocateBorder);
    REGISTER_LOCATOR_BUILDER(SelfLocate2T);
    REGISTER_LOCATOR_BUILDER(SelfLocateLT);
    REGISTER_LOCATOR_BUILDER(SelfLocatePT);
    REGISTER_LOCATOR_BUILDER(SelfLocate2X);
}

/* ------------------------ Rerun log gating (BT) ------------------------ */
// “把所有log日志都套一层 enableLog 判断”
#define RR_ENABLED(BRAIN_PTR) ((BRAIN_PTR) && (BRAIN_PTR)->locator && (BRAIN_PTR)->locator->enableLog && (BRAIN_PTR)->log)

#define RR_LOG(BRAIN_PTR, PATH, OBJ) \
    do { \
        if (RR_ENABLED(BRAIN_PTR)) { \
            (BRAIN_PTR)->log->setTimeNow(); \
            (BRAIN_PTR)->log->log((PATH), (OBJ)); \
        } \
    } while (0)

/* ------------------------ Locator ------------------------ */

void Locator::init(const FieldDimensions& fd,
                   int minMarkerCntParam,
                   double residualToleranceParam,
                   double muOffestParam,
                   bool enableLogParam,
                   string logIPParam)
{
    fieldDimensions = fd;
    calcFieldMarkers(fd); // 内部会 rebuildFieldMarkerIndex()

    minMarkerCnt = minMarkerCntParam;
    residualTolerance = residualToleranceParam;
    muOffset = muOffestParam;
    enableLog = enableLogParam;
    logIP = std::move(logIPParam);

    if (enableLog) {
        auto connectError = log.connect(logIP);
        if (connectError.is_err()) prtErr(format("Rerun log connect Error: %s", connectError.description.c_str()));
        auto saveError = log.save("/home/booster/log.rrd");
        if (saveError.is_err()) prtErr(format("Rerun log save Error: %s", saveError.description.c_str()));
    }
}

void Locator::calcFieldMarkers(const FieldDimensions& fd)
{
    // 避免 init() 多次调用导致 markers 累积 & 迭代变慢
    fieldMarkers.clear();
    fieldMarkers.reserve(26);

    fieldMarkers.emplace_back(FieldMarker{'X', 0.0, -fd.circleRadius, 0.0});
    fieldMarkers.emplace_back(FieldMarker{'X', 0.0,  fd.circleRadius, 0.0});

    fieldMarkers.emplace_back(FieldMarker{'P',  fd.length / 2 - fd.penaltyDist, 0.0, 0.0});
    fieldMarkers.emplace_back(FieldMarker{'P', -fd.length / 2 + fd.penaltyDist, 0.0, 0.0});

    fieldMarkers.emplace_back(FieldMarker{'T', 0.0,  fd.width / 2, 0.0});
    fieldMarkers.emplace_back(FieldMarker{'T', 0.0, -fd.width / 2, 0.0});

    fieldMarkers.emplace_back(FieldMarker{'L',  (fd.length / 2 - fd.penaltyAreaLength),  fd.penaltyAreaWidth / 2, 0.0});
    fieldMarkers.emplace_back(FieldMarker{'L',  (fd.length / 2 - fd.penaltyAreaLength), -fd.penaltyAreaWidth / 2, 0.0});
    fieldMarkers.emplace_back(FieldMarker{'L', -(fd.length / 2 - fd.penaltyAreaLength),  fd.penaltyAreaWidth / 2, 0.0});
    fieldMarkers.emplace_back(FieldMarker{'L', -(fd.length / 2 - fd.penaltyAreaLength), -fd.penaltyAreaWidth / 2, 0.0});
    fieldMarkers.emplace_back(FieldMarker{'T',  fd.length / 2,  fd.penaltyAreaWidth / 2, 0.0});
    fieldMarkers.emplace_back(FieldMarker{'T',  fd.length / 2, -fd.penaltyAreaWidth / 2, 0.0});
    fieldMarkers.emplace_back(FieldMarker{'T', -fd.length / 2,  fd.penaltyAreaWidth / 2, 0.0});
    fieldMarkers.emplace_back(FieldMarker{'T', -fd.length / 2, -fd.penaltyAreaWidth / 2, 0.0});

    fieldMarkers.emplace_back(FieldMarker{'L',  (fd.length / 2 - fd.goalAreaLength),  fd.goalAreaWidth / 2, 0.0});
    fieldMarkers.emplace_back(FieldMarker{'L',  (fd.length / 2 - fd.goalAreaLength), -fd.goalAreaWidth / 2, 0.0});
    fieldMarkers.emplace_back(FieldMarker{'L', -(fd.length / 2 - fd.goalAreaLength),  fd.goalAreaWidth / 2, 0.0});
    fieldMarkers.emplace_back(FieldMarker{'L', -(fd.length / 2 - fd.goalAreaLength), -fd.goalAreaWidth / 2, 0.0});
    fieldMarkers.emplace_back(FieldMarker{'T',  fd.length / 2,  fd.goalAreaWidth / 2, 0.0});
    fieldMarkers.emplace_back(FieldMarker{'T',  fd.length / 2, -fd.goalAreaWidth / 2, 0.0});
    fieldMarkers.emplace_back(FieldMarker{'T', -fd.length / 2,  fd.goalAreaWidth / 2, 0.0});
    fieldMarkers.emplace_back(FieldMarker{'T', -fd.length / 2, -fd.goalAreaWidth / 2, 0.0});

    fieldMarkers.emplace_back(FieldMarker{'L',  fd.length / 2,  fd.width / 2, 0.0});
    fieldMarkers.emplace_back(FieldMarker{'L',  fd.length / 2, -fd.width / 2, 0.0});
    fieldMarkers.emplace_back(FieldMarker{'L', -fd.length / 2,  fd.width / 2, 0.0});
    fieldMarkers.emplace_back(FieldMarker{'L', -fd.length / 2, -fd.width / 2, 0.0});

    // 按 type 建索引：minDist/getOffset 只扫同类点
    rebuildFieldMarkerIndex();
}

vector<Locator::ObsMarker> Locator::buildObsMarkers(const vector<FieldMarker>& markers_r) const
{
    vector<ObsMarker> obs;
    obs.reserve(markers_r.size());

    for (const auto& m : markers_r)
    {
        const double dist = std::max(norm(m.x, m.y), 0.1);
        obs.push_back(ObsMarker{m.type, m.x, m.y, 3.0 / dist});
    }
    return obs;
}

int Locator::genInitialParticles(int num)
{
    // NOTE: 这里的 srand 并不会影响 Eigen::Random()；保留原逻辑不动
    unsigned long long seed =
        chr::duration_cast<std::chrono::milliseconds>(chr::system_clock::now().time_since_epoch()).count();
    srand(static_cast<unsigned int>(seed));

    hypos.resize(num, 6);
    hypos.leftCols<3>() = Eigen::MatrixXd::Random(num, 3);

    auto [xmin, xmax, ymin, ymax, thetamin, thetamax] = constraints;
    hypos.col(0) = hypos.col(0) * (xmax - xmin) / 2.0 + (xmin + xmax) / 2.0;
    hypos.col(1) = hypos.col(1) * (ymax - ymin) / 2.0 + (ymin + ymax) / 2.0;
    hypos.col(2) = hypos.col(2) * (thetamax - thetamin) / 2.0 + (thetamin + thetamax) / 2.0;

    logParticles();
    return 0;
}

int Locator::genParticles()
{
    const int oldN = static_cast<int>(hypos.rows());
    int num = static_cast<int>(oldN * numShrinkRatio);
    if (num <= 0 || oldN <= 0) return 1;

    // 写入 hypos_buf，然后 swap，避免 old_hypos = hypos 的整拷贝
    hypos_buf.resize(num, 6);

    // -----------------------------
    // 重采样（Systematic Resampling，O(N)）
    // 依赖 hypos.col(5) 为 CDF（单调递增，最后≈1）
    // -----------------------------
    const double invN = 1.0 / static_cast<double>(num);
    const double u0 = ((Eigen::VectorXd::Random(1)(0) + 1.0) * 0.5) * invN; // [0, 1/N)

    int j = 0;
    for (int i = 0; i < num; ++i)
    {
        const double u = u0 + static_cast<double>(i) * invN;
        while (j < oldN - 1 && hypos(j, 5) < u) ++j;
        hypos_buf.row(i).head<3>() = hypos.row(j).head<3>();
    }

    // swap：新粒子成为 hypos
    hypos.swap(hypos_buf);

    // 添加扰动（逐步收缩）
    offsetX *= offsetShrinkRatio;
    offsetY *= offsetShrinkRatio;
    offsetTheta *= offsetShrinkRatio;

    Eigen::MatrixXd offsets = Eigen::MatrixXd::Random(num, 3);
    offsets.col(0) *= offsetX;
    offsets.col(1) *= offsetY;
    offsets.col(2) *= offsetTheta;

    hypos.leftCols<3>() += offsets;

    // clip to constraints
    hypos.col(0) = hypos.col(0).cwiseMax(constraints.xmin).cwiseMin(constraints.xmax);
    hypos.col(1) = hypos.col(1).cwiseMax(constraints.ymin).cwiseMin(constraints.ymax);
    hypos.col(2) = hypos.col(2).cwiseMax(constraints.thetamin).cwiseMin(constraints.thetamax);

    logParticles();
    return 0;
}

double Locator::minDistXY(char type, double x, double y) const
{
    const auto& targets = fieldMarkersByType[static_cast<unsigned char>(type)];
    if (targets.empty()) return std::numeric_limits<double>::infinity();

    double minD2 = std::numeric_limits<double>::infinity();
    for (const auto& t : targets)
    {
        const double dx = t.x - x;
        const double dy = t.y - y;
        const double d2 = dx * dx + dy * dy;
        if (d2 < minD2) minD2 = d2;
    }
    return std::sqrt(minD2);
}

array<double, 2> Locator::getOffsetXY(char type, double x, double y) const
{
    const auto& targets = fieldMarkersByType[static_cast<unsigned char>(type)];
    if (targets.empty()) return {0.0, 0.0};

    double minD2 = std::numeric_limits<double>::infinity();
    const FieldMarker* best = nullptr;

    for (const auto& t : targets)
    {
        const double dx = t.x - x;
        const double dy = t.y - y;
        const double d2 = dx * dx + dy * dy;
        if (d2 < minD2) { minD2 = d2; best = &t; }
    }
    if (!best) return {0.0, 0.0};
    return {best->x - x, best->y - y};
}

double Locator::minDist(const FieldMarker& marker)
{
    return minDistXY(marker.type, marker.x, marker.y);
}

array<double, 2> Locator::getOffset(const FieldMarker& marker)
{
    return getOffsetXY(marker.type, marker.x, marker.y);
}

inline double Locator::residual(const vector<ObsMarker>& obs, const Pose2D& pose)
{
    double res = 0.0;

    // tight loop：对同一个 pose，只算一次 sin/cos
    const double c = std::cos(pose.theta);
    const double s = std::sin(pose.theta);

    for (const auto& o : obs)
    {
        // 不构造 FieldMarker，直接算 field frame 坐标
        const double x_f = c * o.x - s * o.y + pose.x;
        const double y_f = s * o.x + c * o.y + pose.y;

        res += minDistXY(o.type, x_f, y_f) * o.inv_dist_w;
    }
    return res;
}

// 保持旧接口（外部验证用）：内部临时 buildObsMarkers
double Locator::residual(const vector<FieldMarker>& markers_r, const Pose2D& pose)
{
    const auto obs = buildObsMarkers(markers_r);
    return residual(obs, pose);
}

Pose2D Locator::finalAdjust(const vector<FieldMarker>& markers_r, const Pose2D& pose)
{
    if (markers_r.empty()) return Pose2D{0, 0, 0};

    double dx = 0.0;
    double dy = 0.0;

    const double c = std::cos(pose.theta);
    const double s = std::sin(pose.theta);

    for (const auto& marker_r : markers_r)
    {
        const double x_f = c * marker_r.x - s * marker_r.y + pose.x;
        const double y_f = s * marker_r.x + c * marker_r.y + pose.y;

        const auto offset = getOffsetXY(marker_r.type, x_f, y_f);
        dx += offset[0];
        dy += offset[1];
    }

    dx /= static_cast<double>(markers_r.size());
    dy /= static_cast<double>(markers_r.size());

    return Pose2D{pose.x + dx, pose.y + dy, pose.theta};
}

int Locator::calcProbs(const vector<FieldMarker>& markers_r)
{
    const auto obs = buildObsMarkers(markers_r);
    return calcProbs(obs);
}

int Locator::calcProbs(const vector<ObsMarker>& obs)
{
    const int rows = static_cast<int>(hypos.rows());
    if (rows <= 1) return 1;

    // -----------------------------
    // 1) 对每个粒子算 residual
    // 2) Welford 单 pass 计算均值/方差
    // 3) 同时跟踪 minRes、bestPose
    // -----------------------------
    double mean = 0.0;
    double M2 = 0.0;
    int n = 0;

    double minRes = std::numeric_limits<double>::infinity();

    for (int i = 0; i < rows; ++i)
    {
        const Pose2D pose{hypos(i, 0), hypos(i, 1), hypos(i, 2)};
        const double res = residual(obs, pose);

        hypos(i, 3) = res;

        // best
        if (res < bestResidual)
        {
            bestResidual = res;
            bestPose = pose;
        }

        // min
        if (res < minRes) minRes = res;

        // Welford update
        ++n;
        const double delta = res - mean;
        mean += delta / static_cast<double>(n);
        const double delta2 = res - mean;
        M2 += delta * delta2;
    }

    const double sigma = (rows > 1) ? std::sqrt(M2 / static_cast<double>(rows - 1)) : 0.0;
    if (sigma < 1e-5) return 1;

    const double mu = minRes - muOffset * sigma;

    // -----------------------------
    // 概率密度常数提取：
    // pdf(r) = (1/(sqrt(2π)*σ)) * exp( -0.5 * ((r-μ)/σ)^2 )
    // -----------------------------
    constexpr double inv_sqrt_2pi = 0.39894228040143267794; // 1/sqrt(2π)
    const double inv_sigma = 1.0 / sigma;
    const double norm_const = inv_sqrt_2pi * inv_sigma;
    const double inv_2sigma2 = 0.5 * inv_sigma * inv_sigma;

    double probSum = 0.0;
    for (int i = 0; i < rows; ++i)
    {
        const double z = (hypos(i, 3) - mu);
        const double p = norm_const * std::exp(-(z * z) * inv_2sigma2);
        hypos(i, 4) = p;
        probSum += p;
    }

    if (std::fabs(probSum) < 1e-12) return 1;

    const double invProbSum = 1.0 / probSum;

    // normalize + build CDF
    double acc = 0.0;
    for (int i = 0; i < rows; ++i)
    {
        const double p = hypos(i, 4) * invProbSum;
        hypos(i, 4) = p;
        acc += p;
        hypos(i, 5) = acc;
    }
    // 防止浮点误差导致末尾 < 1
    hypos(rows - 1, 5) = 1.0;

    return 0;
}

bool Locator::isConverged()
{
    return (
        (hypos.col(0).maxCoeff() - hypos.col(0).minCoeff() < convergeTolerance) &&
        (hypos.col(1).maxCoeff() - hypos.col(1).minCoeff() < convergeTolerance) &&
        (hypos.col(2).maxCoeff() - hypos.col(2).minCoeff() < convergeTolerance)
    );
}

LocateResult Locator::locateRobot(const vector<FieldMarker>& markers_r,
                                 PoseBox2D constraintsParam,
                                 int numParticles,
                                 double offsetXParam,
                                 double offsetYParam,
                                 double offsetThetaParam)
{
    auto start_time = chr::high_resolution_clock::now();
    LocateResult res;

    const double bigEnoughNum = 1e6;
    res.residual = bigEnoughNum;
    bestResidual = bigEnoughNum;
    bestPose = Pose2D{0.0, 0.0, 0.0};

    if (markers_r.size() < static_cast<size_t>(minMarkerCnt))
    {
        res.success = false;
        res.code = 4;
        res.msecs = msecsSince(start_time);
        return res;
    }

    constraints = constraintsParam;
    offsetX = offsetXParam;
    offsetY = offsetYParam;
    offsetTheta = offsetThetaParam;

    // 预计算观测权重（inv_dist_w）：主循环不再重复 sqrt/div
    const auto obs = buildObsMarkers(markers_r);

    genInitialParticles(numParticles);
    if (calcProbs(obs))
    {
        res.success = false;
        res.code = 5;
        res.msecs = msecsSince(start_time);
        return res;
    }

    for (int i = 0; i < maxIteration; i++)
    {
        res.residual = bestResidual / static_cast<double>(markers_r.size());

        if (isConverged())
        {
            if (res.residual > residualTolerance)
            {
                res.success = false;
                res.code = 2;
                res.msecs = msecsSince(start_time);
                return res;
            }

            res.success = true;
            res.code = 0;
            res.pose = bestPose;
            res.pose.theta = toPInPI(res.pose.theta);
            res.msecs = msecsSince(start_time);
            return res;
        }

        if (genParticles())
        {
            res.success = false;
            res.code = 1;
            res.msecs = msecsSince(start_time);
            return res;
        }

        if (calcProbs(obs))
        {
            res.success = false;
            res.code = 5;
            res.msecs = msecsSince(start_time);
            return res;
        }
    }

    res.success = false;
    res.code = 3;
    res.msecs = msecsSince(start_time);
    return res;
}

void Locator::logParticles()
{
    if (!enableLog) return;

    vector<rerun::Position2D> points;
    points.reserve(static_cast<size_t>(hypos.rows()));

    for (int i = 0; i < static_cast<int>(hypos.rows()); i++)
    {
        points.emplace_back(rerun::Position2D{
            static_cast<float>(hypos(i, 0)),
            static_cast<float>(hypos(i, 1))
        });
    }

    log.log(
        "field/hypos",
        rerun::Points2D(points)
            .with_draw_order(20.0)
            .with_colors({rerun::Color{255, 0, 0, 255}})
            .with_radii({0.05})
            .with_draw_order(20)
    );
}

/* ------------------------ Behavior Tree Nodes ------------------------ */

NodeStatus SelfLocate::tick()
{
    auto rr_debug = [&](const string& msg) {
        RR_LOG(brain, "debug/SelfLocate", rerun::TextLog(msg));
    };

    double interval = getInput<double>("msecs_interval").value();
    if (brain->msecsSince(brain->data->lastSuccessfulLocalizeTime) < interval) return NodeStatus::SUCCESS;

    string mode = getInput<string>("mode").value();
    double xMin, xMax, yMin, yMax, thetaMin, thetaMax;
    auto markers = brain->data->getMarkersForLocator();

    if (mode == "face_forward")
    {
        xMin = -brain->config->fieldDimensions.length / 2;
        xMax =  brain->config->fieldDimensions.length / 2;
        yMin = -brain->config->fieldDimensions.width / 2;
        yMax =  brain->config->fieldDimensions.width / 2;
        thetaMin = -M_PI / 4;
        thetaMax =  M_PI / 4;
    }
    else if (mode == "trust_direction")
    {
        int msec = static_cast<int>(brain->msecsSince(brain->data->lastSuccessfulLocalizeTime));
        double maxDriftSpeed = 0.1;
        double maxDrift = msec / 1000.0 * maxDriftSpeed;

        xMin = max(-brain->config->fieldDimensions.length / 2 - 2, brain->data->robotPoseToField.x - maxDrift);
        xMax = min( brain->config->fieldDimensions.length / 2 + 2, brain->data->robotPoseToField.x + maxDrift);
        yMin = max(-brain->config->fieldDimensions.width / 2 - 2, brain->data->robotPoseToField.y - maxDrift);
        yMax = min( brain->config->fieldDimensions.width / 2 + 2, brain->data->robotPoseToField.y + maxDrift);
        thetaMin = brain->data->robotPoseToField.theta - M_PI / 180;
        thetaMax = brain->data->robotPoseToField.theta + M_PI / 180;
    }
    else // fall_recovery
    {
        xMin = -brain->config->fieldDimensions.length / 2 - 2;
        xMax =  brain->config->fieldDimensions.length / 2 + 2;
        yMin = -brain->config->fieldDimensions.width / 2 - 2;
        yMax =  brain->config->fieldDimensions.width / 2 + 2;
        thetaMin = brain->data->robotPoseToField.theta - M_PI / 180;
        thetaMax = brain->data->robotPoseToField.theta + M_PI / 180;
    }

    PoseBox2D constraints{xMin, xMax, yMin, yMax, thetaMin, thetaMax};
    auto res = brain->locator->locateRobot(markers, constraints);

    // 构造 marker 描述字符串：只在需要日志时做（避免无谓开销）
    string mstring;
    if (RR_ENABLED(brain)) {
        mstring.reserve(markers.size() * 32);
        for (int i = 0; i < static_cast<int>(markers.size()); i++) {
            const auto& m = markers[i];
            mstring += format("type: %c  x: %.1f y: %.1f", m.type, m.x, m.y);
        }
    }

    if (res.success) {
        RR_LOG(
            brain,
            "field/recal",
            rerun::Arrows2D::from_vectors({{res.pose.x - brain->data->robotPoseToField.x, -res.pose.y + brain->data->robotPoseToField.y}})
                .with_origins({{brain->data->robotPoseToField.x, -brain->data->robotPoseToField.y}})
                .with_colors(0x00FF00FF)
                .with_radii(0.01)
                .with_draw_order(10)
                .with_labels({"pf"})
        );
    }

    rr_debug(
        format(
            "success: %d  residual: %.2f  marker.size: %d  minMarkerCnt: %d  resTolerance: %.2f marker: %s",
            res.success,
            res.residual,
            static_cast<int>(markers.size()),
            brain->locator->minMarkerCnt,
            brain->locator->residualTolerance,
            mstring.c_str()
        )
    );

    if (!res.success) return NodeStatus::SUCCESS;

    brain->calibrateOdom(res.pose.x, res.pose.y, res.pose.theta);
    brain->tree->setEntry<bool>("odom_calibrated", true);
    brain->data->lastSuccessfulLocalizeTime = brain->get_clock()->now();
    prtDebug("定位成功: " + to_string(res.pose.x) + " " + to_string(res.pose.y) + " " + to_string(rad2deg(res.pose.theta)) + " Dur: " + to_string(res.msecs));

    return NodeStatus::SUCCESS;
}

NodeStatus SelfLocateEnterField::tick()
{
    auto rr_log = [&](const string& msg, bool success) {
        RR_LOG(
            brain,
            "debug/SelfLocateEnterField",
            rerun::TextLog(msg).with_level(success ? rerun::TextLogLevel::Info : rerun::TextLogLevel::Error)
        );
    };

    double interval = getInput<double>("msecs_interval").value();
    if (brain->msecsSince(brain->data->lastSuccessfulLocalizeTime) < interval) return NodeStatus::SUCCESS;

    auto markers = brain->data->getMarkersForLocator();
    auto fd = brain->config->fieldDimensions;

    PoseBox2D cEnterLeft  = {-fd.length / 2, -fd.circleRadius,  fd.width / 2,  fd.width / 2 + 1, -M_PI / 2 - M_PI / 6, -M_PI / 2 + M_PI / 6};
    PoseBox2D cEnterRight = {-fd.length / 2, -fd.circleRadius, -fd.width / 2 - 1, -fd.width / 2,      M_PI / 2 - M_PI / 6,  M_PI / 2 + M_PI / 6};

    auto resLeft = brain->locator->locateRobot(markers, cEnterLeft);
    auto resRight = brain->locator->locateRobot(markers, cEnterRight);
    LocateResult res;

    static string lastReport = "";
    string report = lastReport;
    if (resLeft.success && !resRight.success) {
        res = resLeft;
        report = "Entering Left";
    }
    else if (!resLeft.success && resRight.success) {
        res = resRight;
        report = "Entering Right";
    }
    else if (resLeft.success && resRight.success) {
        if (resLeft.residual < resRight.residual) {
            res = resLeft;
            report = "Entering Left";
        }
        else {
            res = resRight;
            report = "Entering Right";
        }
    } else {
        res = resLeft;
    }

    if (report != lastReport) {
        brain->speak(report);
        lastReport = report;
    }

    rr_log(
        format(
            "%s left success: %d  left residual: %.2f  right success %d  right residual %.2f resTolerance: %.2f markers: %d minMarkerCnt: %d ",
            report.c_str(),
            resLeft.success,
            resLeft.residual,
            resRight.success,
            resRight.residual,
            brain->locator->residualTolerance,
            static_cast<int>(markers.size()),
            brain->locator->minMarkerCnt
        ),
        res.success
    );

    RR_LOG(
        brain,
        "field/recal_enter_field",
        rerun::Arrows2D::from_vectors({{res.pose.x - brain->data->robotPoseToField.x, -res.pose.y + brain->data->robotPoseToField.y}})
            .with_origins({{brain->data->robotPoseToField.x, -brain->data->robotPoseToField.y}})
            .with_colors(res.success ? 0x00FF00FF : 0xFF0000FF)
            .with_radii(0.01)
            .with_draw_order(10)
            .with_labels({"pfe"})
    );

    if (!res.success) return NodeStatus::SUCCESS;

    brain->calibrateOdom(res.pose.x, res.pose.y, res.pose.theta);
    brain->tree->setEntry<bool>("odom_calibrated", true);
    brain->data->lastSuccessfulLocalizeTime = brain->get_clock()->now();
    prtDebug("定位成功: " + to_string(res.pose.x) + " " + to_string(res.pose.y) + " " + to_string(rad2deg(res.pose.theta)) + " Dur: " + to_string(res.msecs));

    return NodeStatus::SUCCESS;
}

NodeStatus SelfLocate1M::tick()
{
    double interval = getInput<double>("msecs_interval").value();
    double maxDist = getInput<double>("max_dist").value();
    if (brain->client->isStandingStill(2000)) maxDist *= 1.5;
    double maxDrift = getInput<double>("max_drift").value();
    bool validate = getInput<bool>("validate").value();

    const string logPathS = "/locate/1m/success";
    const string logPathF = "/locate/1m/fail";

    auto msecs = brain->msecsSince(brain->data->lastSuccessfulLocalizeTime);
    if (msecs < interval) {
        RR_LOG(brain, logPathF, rerun::TextLog(format("Failed, msecs(%.1f) < interval(%.1f)", msecs, interval)));
        return NodeStatus::SUCCESS;
    }

    int markerIndex = -1;
    GameObject marker;
    MapMarking mapMarker;
    double minDist = 100.0;
    auto markings = brain->data->getMarkings();
    for (int i = 0; i < static_cast<int>(markings.size()); i++) {
        const auto& m = markings[i];
        if (m.name == "LOLG" || m.name == "LORG" || m.name == "LSLG" || m.name == "LSRG") continue;

        if (m.range < minDist) {
            minDist = m.range;
            markerIndex = i;
            marker = m;
        }
    }

    if (markerIndex < 0 || markerIndex >= static_cast<int>(markings.size())
        || marker.id < 0 || marker.id >= static_cast<int>(brain->config->mapMarkings.size()))
    {
        RR_LOG(brain, logPathF, rerun::TextLog("Failed, No markings Found. Or marker id invalid."));
        return NodeStatus::SUCCESS;
    }
    mapMarker = brain->config->mapMarkings[marker.id];

    if (marker.range > maxDist) {
        RR_LOG(brain, logPathF, rerun::TextLog(format("Failed, min marker Dist(%.2f) > maxDist(%.2f)", marker.range, maxDist)));
        return NodeStatus::SUCCESS;
    }

    if (!brain->isBoundingBoxInCenter(marker.boundingBox)) {
        RR_LOG(brain, logPathF, rerun::TextLog("Failed, boundingbox is not in the center area"));
        return NodeStatus::SUCCESS;
    }

    double dx = mapMarker.x - marker.posToField.x;
    double dy = mapMarker.y - marker.posToField.y;

    double drift = norm(dx, dy);
    if (drift > maxDrift) {
        RR_LOG(brain, logPathF, rerun::TextLog(format("Failed, drift(%.2f) > maxDrift(%.2f)", drift, maxDrift)));
        return NodeStatus::SUCCESS;
    }

    Pose2D hypoPose = brain->data->robotPoseToField;
    hypoPose.x += dx;
    hypoPose.y += dy;

    auto allMarkers = brain->data->getMarkersForLocator();
    if (validate && !allMarkers.empty()) {
        double residual = brain->locator->residual(allMarkers, hypoPose) / allMarkers.size();
        if (residual > brain->locator->residualTolerance) {
            RR_LOG(brain, logPathF, rerun::TextLog(format("Failed, validation residual(%.2f) > tolerance(%.2f)", residual, brain->locator->residualTolerance)));
            return NodeStatus::SUCCESS;
        }
    }

    RR_LOG(brain, logPathS, rerun::TextLog(format("Success. Drift = %.2f", drift)));
    RR_LOG(
        brain,
        "field/recal/1m/success",
        rerun::Arrows2D::from_vectors({{hypoPose.x - brain->data->robotPoseToField.x, -hypoPose.y + brain->data->robotPoseToField.y}})
            .with_origins({{brain->data->robotPoseToField.x, -brain->data->robotPoseToField.y}})
            .with_colors(0x00FF00FF)
            .with_radii(0.01)
            .with_draw_order(10)
            .with_labels({marker.name})
    );

    brain->calibrateOdom(hypoPose.x, hypoPose.y, hypoPose.theta);
    brain->data->lastSuccessfulLocalizeTime = brain->get_clock()->now();
    return NodeStatus::SUCCESS;
}

NodeStatus SelfLocate2X::tick()
{
    double interval = getInput<double>("msecs_interval").value();
    double maxDist = getInput<double>("max_dist").value();
    if (brain->client->isStandingStill(2000)) maxDist *= 1.5;
    double maxDrift = getInput<double>("max_drift").value();
    bool validate = getInput<bool>("validate").value();

    const string logPathS = "/locate/2x/success";
    const string logPathF = "/locate/2x/fail";

    auto msecs = brain->msecsSince(brain->data->lastSuccessfulLocalizeTime);
    if (msecs < interval) {
        RR_LOG(brain, logPathF, rerun::TextLog(format("Failed, msecs(%.1f) < interval(%.1f)", msecs, interval)));
        return NodeStatus::SUCCESS;
    }

    auto points = brain->data->getMarkingsByType({"XCross"});
    if (points.size() != 2) {
        RR_LOG(brain, logPathF, rerun::TextLog(format("Failed, point cnt(%d) != 2", static_cast<int>(points.size()))));
        return NodeStatus::SUCCESS;
    }

    const auto& p0 = points[0];
    const auto& p1 = points[1];

    if (p0.range > maxDist || p1.range > maxDist) {
        RR_LOG(brain, logPathF, rerun::TextLog(format("Failed, p0 range (%.2f) or p1 range (%.2f) > maxDist(%.2f)", p0.range, p1.range, maxDist)));
        return NodeStatus::SUCCESS;
    }

    double xDist = fabs(p0.posToField.x - p1.posToField.x);
    if (xDist > 0.5) {
        RR_LOG(brain, logPathF, rerun::TextLog(format("Failed, xDist(%.2f) > maxDist(%.2f)", xDist, 0.5)));
        return NodeStatus::SUCCESS;
    }

    double yDist = fabs(p0.posToField.y - p1.posToField.y);
    double mapYDist = brain->config->fieldDimensions.circleRadius * 2.0;
    if (fabs(yDist - mapYDist) > 0.5) {
        RR_LOG(brain, logPathF, rerun::TextLog(format("Failed, yDist(%.2f) too far (%.2f) from mapYDist(%.2f)", yDist, 0.5, mapYDist)));
        return NodeStatus::SUCCESS;
    }

    double dx = - (p0.posToField.x + p1.posToField.x) / 2.0;
    double dy = - (p0.posToField.y + p1.posToField.y) / 2.0;
    double drift = norm(dx, dy);

    if (drift > maxDrift) {
        RR_LOG(brain, logPathF, rerun::TextLog(format("Failed, dirft(%.2f) > maxDrift(%.2f)", drift, maxDrift)));
        return NodeStatus::SUCCESS;
    }

    Pose2D hypoPose = brain->data->robotPoseToField;
    hypoPose.x += dx;
    hypoPose.y += dy;

    auto allMarkers = brain->data->getMarkersForLocator();
    if (validate && !allMarkers.empty()) {
        double residual = brain->locator->residual(allMarkers, hypoPose) / allMarkers.size();
        if (residual > brain->locator->residualTolerance) {
            RR_LOG(brain, logPathF, rerun::TextLog(format("Failed, validation residual(%.2f) > tolerance(%.2f)", residual, brain->locator->residualTolerance)));
            return NodeStatus::SUCCESS;
        }
    }

    RR_LOG(brain, logPathS, rerun::TextLog(format("Success. Dist = %.2f", drift)));
    RR_LOG(
        brain,
        "field/recal/2x/success",
        rerun::Arrows2D::from_vectors({{hypoPose.x - brain->data->robotPoseToField.x, -hypoPose.y + brain->data->robotPoseToField.y}})
            .with_origins({{brain->data->robotPoseToField.x, -brain->data->robotPoseToField.y}})
            .with_colors(0x00FF00FF)
            .with_radii(0.01)
            .with_draw_order(10)
            .with_labels({"1p"})
    );

    brain->calibrateOdom(hypoPose.x, hypoPose.y, hypoPose.theta);
    brain->data->lastSuccessfulLocalizeTime = brain->get_clock()->now();
    return NodeStatus::SUCCESS;
}

NodeStatus SelfLocate2T::tick()
{
    double interval = getInput<double>("msecs_interval").value();
    double maxDist = getInput<double>("max_dist").value();
    if (brain->client->isStandingStill(2000)) maxDist *= 1.5;
    double maxDrift = getInput<double>("max_drift").value();
    bool validate = getInput<bool>("validate").value();

    const string logPathS = "/locate/2t/success";
    const string logPathF = "/locate/2t/fail";

    auto msecs = brain->msecsSince(brain->data->lastSuccessfulLocalizeTime);
    if (msecs < interval) {
        RR_LOG(brain, logPathF, rerun::TextLog(format("Failed, msecs(%.1f) < interval(%.1f)", msecs, interval)));
        return NodeStatus::SUCCESS;
    }

    auto markers = brain->data->getMarkingsByType({"TCross"});
    GameObject m1, m2;
    bool found = false;
    auto fd = brain->config->fieldDimensions;

    for (int i = 0; i < static_cast<int>(markers.size()); i++) {
        m1 = markers[i];
        if (m1.range > maxDist) continue;

        for (int j = i + 1; j < static_cast<int>(markers.size()); j++) {
            m2 = markers[j];
            if (m2.range > maxDist) continue;

            if (fabs(m1.posToField.x - m2.posToField.x) < 0.3
                && fabs(fabs(m1.posToField.y - m2.posToField.y) - fabs(fd.goalAreaWidth - fd.penaltyAreaWidth) / 2.0) < 0.3)
            {
                found = true;
                break;
            }
        }
        if (found) break;
    }

    if (!found) {
        RR_LOG(brain, logPathF, rerun::TextLog(format("Failed, No pattern within maxDist(%.2f) Found", maxDist)));
        return NodeStatus::SUCCESS;
    }

    Point2D pos_o = {
        (m1.posToField.x + m2.posToField.x) / 2,
        (m1.posToField.y + m2.posToField.y) / 2
    };
    Point2D pos_m;

    static constexpr array<double,2> halfs{-1.0, 1.0};
    static constexpr array<double,2> sides{-1.0, 1.0};

    bool matched = false;
    for (auto half : halfs) {
        for (auto side : sides) {
            pos_m = {
                half * (fd.length / 2.0),
                side * (fd.penaltyAreaWidth + fd.goalAreaWidth) / 4.0
            };
            double dist = norm(pos_o.x - pos_m.x, pos_o.y - pos_m.y);
            if (dist < maxDrift) { matched = true; break; }
        }
        if (matched) break;
    }

    if (!matched) {
        RR_LOG(brain, logPathF, rerun::TextLog(format("Failed, can not match to any map positions within maxDrift(%.2f)", maxDrift)));
        return NodeStatus::SUCCESS;
    }

    double dx = pos_m.x - pos_o.x;
    double dy = pos_m.y - pos_o.y;
    double drift = norm(dx, dy);

    Pose2D hypoPose = brain->data->robotPoseToField;
    hypoPose.x += dx;
    hypoPose.y += dy;

    auto allMarkers = brain->data->getMarkersForLocator();
    if (validate && !allMarkers.empty()) {
        double residual = brain->locator->residual(allMarkers, hypoPose) / allMarkers.size();
        if (residual > brain->locator->residualTolerance) {
            RR_LOG(brain, logPathF, rerun::TextLog(format("Failed, validation residual(%.2f) > tolerance(%.2f)", residual, brain->locator->residualTolerance)));
            return NodeStatus::SUCCESS;
        }
    }

    RR_LOG(brain, logPathS, rerun::TextLog(format("Success. Dist = %.2f", drift)));
    RR_LOG(
        brain,
        "field/recal/2t/success",
        rerun::Arrows2D::from_vectors({{hypoPose.x - brain->data->robotPoseToField.x, -hypoPose.y + brain->data->robotPoseToField.y}})
            .with_origins({{brain->data->robotPoseToField.x, -brain->data->robotPoseToField.y}})
            .with_colors(0x00FF00FF)
            .with_radii(0.01)
            .with_draw_order(10)
            .with_labels({"2t"})
    );

    brain->calibrateOdom(hypoPose.x, hypoPose.y, hypoPose.theta);
    brain->data->lastSuccessfulLocalizeTime = brain->get_clock()->now();
    return NodeStatus::SUCCESS;
}

NodeStatus SelfLocateLT::tick()
{
    double interval = getInput<double>("msecs_interval").value();
    double maxDist = getInput<double>("max_dist").value();
    if (brain->client->isStandingStill(2000)) maxDist *= 1.5;
    double maxDrift = getInput<double>("max_drift").value();
    bool validate = getInput<bool>("validate").value();

    const string logPathS = "/locate/lt/success";
    const string logPathF = "/locate/lt/fail";

    auto msecs = brain->msecsSince(brain->data->lastSuccessfulLocalizeTime);
    if (msecs < interval) {
        RR_LOG(brain, logPathF, rerun::TextLog(format("Failed, msecs(%.1f) < interval(%.1f)", msecs, interval)));
        return NodeStatus::SUCCESS;
    }

    auto tMarkers = brain->data->getMarkingsByType({"TCross"});
    auto lMarkers = brain->data->getMarkingsByType({"LCross"});

    GameObject t, l;
    bool found = false;
    auto fd = brain->config->fieldDimensions;

    for (int i = 0; i < static_cast<int>(tMarkers.size()); i++) {
        t = tMarkers[i];
        if (t.range > maxDist) continue;

        for (int j = 0; j < static_cast<int>(lMarkers.size()); j++) {
            l = lMarkers[j];
            if (l.range > maxDist) continue;

            if (fabs(t.posToField.y - l.posToField.y) < 0.3
                && fabs(fabs(t.posToField.x - l.posToField.x) - fd.goalAreaLength) < 0.3)
            {
                found = true;
                break;
            }
        }
        if (found) break;
    }

    if (!found) {
        RR_LOG(brain, logPathF, rerun::TextLog(format("Failed, No pattern within MaxDist(%.2f) Found", maxDist)));
        return NodeStatus::SUCCESS;
    }

    Point2D pos_o = {
        (t.posToField.x + l.posToField.x) / 2,
        (t.posToField.y + l.posToField.y) / 2
    };
    Point2D pos_m;

    static constexpr array<double,2> halfs{-1.0, 1.0};
    static constexpr array<double,2> sides{-1.0, 1.0};

    bool matched = false;
    for (auto half : halfs) {
        for (auto side : sides) {
            pos_m = {
                half * (fd.length / 2.0 - fd.goalAreaLength / 2.0),
                side * (fd.goalAreaWidth / 2.0)
            };
            double dist = norm(pos_o.x - pos_m.x, pos_o.y - pos_m.y);
            if (dist < maxDrift) { matched = true; break; }
        }
        if (matched) break;
    }

    if (!matched) {
        RR_LOG(brain, logPathF, rerun::TextLog(format("Failed, can not match to any map positions within maxDrift(%.2f)", maxDrift)));
        return NodeStatus::SUCCESS;
    }

    double dx = pos_m.x - pos_o.x;
    double dy = pos_m.y - pos_o.y;
    double drift = norm(dx, dy);

    Pose2D hypoPose = brain->data->robotPoseToField;
    hypoPose.x += dx;
    hypoPose.y += dy;

    auto allMarkers = brain->data->getMarkersForLocator();
    if (validate && !allMarkers.empty()) {
        double residual = brain->locator->residual(allMarkers, hypoPose) / allMarkers.size();
        if (residual > brain->locator->residualTolerance) {
            RR_LOG(brain, logPathF, rerun::TextLog(format("Failed, validation residual(%.2f) > tolerance(%.2f)", residual, brain->locator->residualTolerance)));
            return NodeStatus::SUCCESS;
        }
    }

    RR_LOG(brain, logPathS, rerun::TextLog(format("Success. Dist = %.2f", drift)));
    RR_LOG(
        brain,
        "field/recal/lt/success",
        rerun::Arrows2D::from_vectors({{hypoPose.x - brain->data->robotPoseToField.x, -hypoPose.y + brain->data->robotPoseToField.y}})
            .with_origins({{brain->data->robotPoseToField.x, -brain->data->robotPoseToField.y}})
            .with_colors(0x00FF00FF)
            .with_radii(0.01)
            .with_draw_order(10)
            .with_labels({"lt"})
    );

    brain->calibrateOdom(hypoPose.x, hypoPose.y, hypoPose.theta);
    brain->data->lastSuccessfulLocalizeTime = brain->get_clock()->now();
    return NodeStatus::SUCCESS;
}

NodeStatus SelfLocatePT::tick()
{
    double interval = getInput<double>("msecs_interval").value();
    double maxDist = getInput<double>("max_dist").value();
    if (brain->client->isStandingStill(2000)) maxDist *= 1.5;
    double maxDrift = getInput<double>("max_drift").value();
    bool validate = getInput<bool>("validate").value();

    const string logPathS = "/locate/pt/success";
    const string logPathF = "/locate/pt/fail";

    auto msecs = brain->msecsSince(brain->data->lastSuccessfulLocalizeTime);
    if (msecs < interval) {
        RR_LOG(brain, logPathF, rerun::TextLog(format("Failed, msecs(%.1f) < interval(%.1f)", msecs, interval)));
        return NodeStatus::SUCCESS;
    }

    auto posts = brain->data->getGoalposts();
    auto tMarkers = brain->data->getMarkingsByType({"TCross"});

    GameObject p, t;
    bool found = false;
    auto fd = brain->config->fieldDimensions;

    for (int i = 0; i < static_cast<int>(posts.size()); i++) {
        p = posts[i];
        if (p.range > maxDist) continue;

        for (int j = 0; j < static_cast<int>(tMarkers.size()); j++) {
            t = tMarkers[j];
            if (t.range > maxDist) continue;

            if (fabs(t.posToField.x - p.posToField.x) < 0.5
                && fabs(fabs(t.posToField.x - p.posToField.x) - fabs(fd.goalAreaWidth - fd.goalWidth) / 2.0) < 0.3)
            {
                found = true;
                break;
            }
        }
        if (found) break;
    }

    if (!found) {
        RR_LOG(brain, logPathF, rerun::TextLog(format("Failed, No pattern within maxDist(%.2f) Found", maxDist)));
        return NodeStatus::SUCCESS;
    }

    Point2D pos_o = {t.posToField.x, t.posToField.y};
    Point2D pos_m;

    static constexpr array<double,2> halfs{-1.0, 1.0};
    static constexpr array<double,2> sides{-1.0, 1.0};

    bool matched = false;
    for (auto half : halfs) {
        for (auto side : sides) {
            pos_m = {half * (fd.length), side * (fd.goalAreaWidth / 2.0)};
            double dist = norm(pos_o.x - pos_m.x, pos_o.y - pos_m.y);
            if (dist < maxDrift) { matched = true; break; }
        }
        if (matched) break;
    }

    if (!matched) {
        RR_LOG(brain, logPathF, rerun::TextLog(format("Failed, can not match to any map positions within maxDrift(%.2f)", maxDrift)));
        return NodeStatus::SUCCESS;
    }

    double dx = pos_m.x - pos_o.x;
    double dy = pos_m.y - pos_o.y;
    double drift = norm(dx, dy);

    Pose2D hypoPose = brain->data->robotPoseToField;
    hypoPose.x += dx;
    hypoPose.y += dy;

    auto allMarkers = brain->data->getMarkersForLocator();
    if (validate && !allMarkers.empty()) {
        double residual = brain->locator->residual(allMarkers, hypoPose) / allMarkers.size();
        if (residual > brain->locator->residualTolerance) {
            RR_LOG(brain, logPathF, rerun::TextLog(format("Failed, validation residual(%.2f) > tolerance(%.2f)", residual, brain->locator->residualTolerance)));
            return NodeStatus::SUCCESS;
        }
    }

    RR_LOG(brain, logPathS, rerun::TextLog(format("Success. Dist = %.2f", drift)));
    RR_LOG(
        brain,
        "field/recal/pt/success",
        rerun::Arrows2D::from_vectors({{hypoPose.x - brain->data->robotPoseToField.x, -hypoPose.y + brain->data->robotPoseToField.y}})
            .with_origins({{brain->data->robotPoseToField.x, -brain->data->robotPoseToField.y}})
            .with_colors(0x00FF00FF)
            .with_radii(0.01)
            .with_draw_order(10)
            .with_labels({"pt"})
    );

    brain->calibrateOdom(hypoPose.x, hypoPose.y, hypoPose.theta);
    brain->data->lastSuccessfulLocalizeTime = brain->get_clock()->now();
    return NodeStatus::SUCCESS;
}

NodeStatus SelfLocateBorder::tick()
{
    double interval = getInput<double>("msecs_interval").value();
    double maxDist = getInput<double>("max_dist").value();
    if (brain->client->isStandingStill(2000)) maxDist *= 1.5;
    double maxDrift = getInput<double>("max_drift").value();
    bool validate = getInput<bool>("validate").value();

    const string logPathS = "/locate/border/success";
    const string logPathF = "/locate/border/fail";

    auto msecs = brain->msecsSince(brain->data->lastSuccessfulLocalizeTime);
    if (msecs < interval) {
        RR_LOG(brain, logPathF, rerun::TextLog(format("Failed, msecs(%.1f) < interval(%.1f)", msecs, interval)));
        return NodeStatus::SUCCESS;
    }

    bool touchLineFound = false;
    FieldLine touchLine;
    double bestConfidenceTouchline = 0.0;

    bool goalLineFound = false;
    FieldLine goalLine;
    double bestConfidenceGoalline = 0.0;

    bool middleLineFound = false;
    FieldLine middleLine;
    double bestConfidenceMiddleLine = 0.0;

    auto fieldLines = brain->data->getFieldLines();
    for (int i = 0; i < static_cast<int>(fieldLines.size()); i++) {
        const auto& line = fieldLines[i];
        if (line.type != LineType::TouchLine && line.type != LineType::GoalLine && line.type != LineType::MiddleLine) continue;
        if (line.confidence < 0.8) continue;

        double dist = pointMinDistToLine(
            Point2D({brain->data->robotPoseToField.x, brain->data->robotPoseToField.y}),
            line.posToField
        );
        if (dist > maxDist) continue;

        if (line.type == LineType::TouchLine) {
            if (line.confidence > bestConfidenceTouchline) {
                bestConfidenceTouchline = line.confidence;
                touchLine = line;
                touchLineFound = true;
            }
        } else if (line.type == LineType::GoalLine) {
            if (line.confidence > bestConfidenceGoalline) {
                bestConfidenceGoalline = line.confidence;
                goalLine = line;
                goalLineFound = true;
            }
        } else if (line.type == LineType::MiddleLine) {
            if (line.confidence > bestConfidenceMiddleLine) {
                bestConfidenceMiddleLine = line.confidence;
                middleLine = line;
                middleLineFound = true;
            }
        }
    }

    double dx = 0.0;
    double dy = 0.0;
    auto fd = brain->config->fieldDimensions;

    if (touchLineFound) {
        double y_m = touchLine.side == LineSide::Left ? fd.width / 2.0 : -fd.width / 2.0;
        double perpDist = pointPerpDistToLine(
            Point2D({brain->data->robotPoseToField.x, brain->data->robotPoseToField.y}),
            touchLine.posToField
        );
        double y_o = touchLine.side == LineSide::Left ?
            brain->data->robotPoseToField.y - perpDist :
            brain->data->robotPoseToField.y + perpDist;
        dy = y_m - y_o;
    }

    if (goalLineFound) {
        double x_m = goalLine.half == LineHalf::Opponent ? fd.length / 2.0 : -fd.length / 2.0;
        double perpDist = pointPerpDistToLine(
            Point2D({brain->data->robotPoseToField.x, brain->data->robotPoseToField.y}),
            goalLine.posToField
        );
        double x_o = goalLine.half == LineHalf::Opponent ?
            brain->data->robotPoseToField.x - perpDist :
            brain->data->robotPoseToField.x + perpDist;
        dx = x_m - x_o;
    }
    else if (middleLineFound) {
        double x_m = 0.0;
        auto linePos = middleLine.posToField;
        auto robotPose = brain->data->robotPoseToField;

        vector<double> pointA(2);
        vector<double> pointB(2);
        vector<double> pointR = {robotPose.x, robotPose.y};

        if (linePos.y0 > linePos.y1) {
            pointA = {linePos.x0, linePos.y0};
            pointB = {linePos.x1, linePos.y1};
        } else {
            pointA = {linePos.x1, linePos.y1};
            pointB = {linePos.x0, linePos.y0};
        }

        vector<double> vl = {pointB[0] - pointA[0], pointB[1] - pointA[1]};
        vector<double> vr = {pointR[0] - pointA[0], pointR[1] - pointA[1]};

        double normvl = norm(vl);
        double normvr = norm(vr);
        if (normvl < 1e-3 || normvr < 1e-3) {
            dx = 10000;
        } else {
            double dist = crossProduct(vr, vl) / normvl;
            double x_o = robotPose.x + dist;
            dx = x_m - x_o;
        }
    }

    if (!touchLineFound && !goalLineFound && !middleLineFound) {
        RR_LOG(brain, logPathF, rerun::TextLog("No touchline or goalline or middleLine found."));
        return NodeStatus::SUCCESS;
    }

    double drift = norm(dx, dy);
    if (drift > maxDrift) {
        RR_LOG(brain, logPathF, rerun::TextLog(format("Failed, drift(%.2f) > maxDrift(%.2f)", drift, maxDrift)));
        return NodeStatus::SUCCESS;
    }

    Pose2D hypoPose = brain->data->robotPoseToField;
    hypoPose.x += dx;
    hypoPose.y += dy;

    auto allMarkers = brain->data->getMarkersForLocator();
    if (validate && !allMarkers.empty()) {
        double residual = brain->locator->residual(allMarkers, hypoPose) / allMarkers.size();
        if (residual > brain->locator->residualTolerance) {
            RR_LOG(brain, logPathF, rerun::TextLog(format("Failed, validation residual(%.2f) > tolerance(%.2f)", residual, brain->locator->residualTolerance)));
            return NodeStatus::SUCCESS;
        }
    }

    RR_LOG(brain, logPathS, rerun::TextLog(format("Success. Drift = %.2f", drift)));

    string label;
    if (touchLineFound) label += "TouchLine";
    if (touchLineFound && (goalLineFound || middleLineFound)) label += " ";
    if (goalLineFound) label += "GoalLine";
    if (middleLineFound) label += "MiddleLine";

    RR_LOG(
        brain,
        "field/recal/border/success",
        rerun::Arrows2D::from_vectors({{hypoPose.x - brain->data->robotPoseToField.x, -hypoPose.y + brain->data->robotPoseToField.y}})
            .with_origins({{brain->data->robotPoseToField.x, -brain->data->robotPoseToField.y}})
            .with_colors(0x00FF00FF)
            .with_radii(0.01)
            .with_draw_order(10)
            .with_labels({label})
    );

    brain->calibrateOdom(hypoPose.x, hypoPose.y, hypoPose.theta);
    brain->data->lastSuccessfulLocalizeTime = brain->get_clock()->now();
    return NodeStatus::SUCCESS;
}